//! Matched consumer for bench.goose. No upstream source changes.
use seqair::bam::{
    CustomizeRecordStore, FilterRawFields, Indel, IndexedBamReader, PileupEngine, Pos0,
    RawBamRecord, RecordStore, record_store::SlimRecord, stream::BamStreamReader,
};
use seqair::reader::IndexPath;
use std::{fs::File, hint::black_box, path::Path, time::Instant};

#[derive(Clone)]
struct Filter;
impl CustomizeRecordStore for Filter {
    type Extra = ();
    fn filter_raw(&mut self, fields: &FilterRawFields<'_>) -> bool {
        fields.flags.raw() & 0x704 == 0
    }
    fn compute(&mut self, _: &SlimRecord, _: &RecordStore) {}
}

fn walk<const DUMP: bool>(engine: &mut PileupEngine) -> (i64, i64) {
    let (mut count, mut sum) = (0, 0);
    while let Some(col) = engine.pileups() {
        let pos = col.pos().as_i64();
        for a in col.alignments() {
            let qpos = a.qpos().map_or(-1, |q| q.as_usize() as i64);
            let op = if a.is_refskip() {
                2
            } else if a.is_del() {
                1
            } else {
                0
            };
            let base = a.base().map_or(b'N', |b| b.as_char() as u8);
            let qual = a.qual().map_or(255, |q| q.as_byte());
            let indel = match a.indel_after() {
                Indel::None => 0,
                Indel::Insertion(n) => i64::from(n),
                Indel::Deletion(n) => -i64::from(n),
            };
            count += 1;
            sum += pos
                + i64::from(a.record_idx().get())
                + qpos
                + i64::from(base)
                + i64::from(qual)
                + indel
                + op;
            if DUMP {
                println!(
                    "{pos}\t{}\t{qpos}\t{op}\t{indel}\t{base}\t{qual}",
                    std::str::from_utf8(a.qname()).unwrap()
                );
            }
        }
    }
    (count, sum)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    if args.get(1).is_some_and(|s| s == "read-sam") {
        let iterations: usize = args[3].parse()?;
        let mut raw = RawBamRecord::default();
        let mut elapsed = 0;
        let (mut records, mut bases) = (0, 0);
        for iteration in 0..=iterations {
            let start = Instant::now();
            let mut reader = seqair::sam::SamStreamReader::new(std::io::BufReader::new(File::open(&args[2])?))?;
            records = 0; bases = 0;
            while reader.read_record(&mut raw)? { records += 1; bases += raw.seq_len(); }
            if iteration > 0 { elapsed += start.elapsed().as_nanos(); }
        }
        println!("{elapsed} {records} {bases}"); return Ok(());
    }
    if args.get(1).is_some_and(|s| s == "write" || s == "write-sam") {
        return write_benchmark(&args);
    }
    assert!(
        (6..=8).contains(&args.len()),
        "usage: seqair-comparison input.bam contig begin end iterations [indexed|dump|mates|indexed-dump|indexed-mates] [index]"
    );
    let path = Path::new(&args[1]);
    let begin: u32 = args[3].parse()?;
    let end: u32 = args[4].parse()?;
    let iterations: usize = args[5].parse()?;
    assert!(begin < end && iterations > 0);
    let span = (Pos0::new(begin).unwrap()..=Pos0::new(end - 1).unwrap()).into();
    let mode = args.get(6).map_or("sequential", String::as_str);
    let indexed = matches!(mode, "indexed" | "indexed-dump" | "indexed-mates");
    let dump = matches!(mode, "dump" | "indexed-dump");
    let mates = matches!(mode, "mates" | "indexed-mates");
    let index = args.get(7).map(IndexPath::new);
    let mut store = RecordStore::new();
    let mut raw = RawBamRecord::default();
    let mut elapsed = [0u128; 3];
    let mut expected = (0, 0, 0, 0);
    for iteration in 0..=iterations {
        let start = Instant::now();
        let (mut records, mut sum) = (0i64, 0i64);
        if !dump && !mates {
            let mut reader = BamStreamReader::new(File::open(path)?)?;
            while reader.read_record(&mut raw)? {
                records += 1;
                sum += raw.seq_len() as i64;
            }
        }
        let read_ns = start.elapsed().as_nanos();
        let start = Instant::now();
        store.clear();
        if indexed {
            let mut reader = match &index {
                Some(index) => IndexedBamReader::open_with_index(path, index)?,
                None => IndexedBamReader::open(path)?,
            };
            let tid = reader.header().tid(&args[2]).ok_or("unknown contig")?;
            reader.fetch_into_customized(tid, span, &mut store, &mut Filter)?;
        } else {
            let mut reader = BamStreamReader::new(File::open(path)?)?;
            let tid = reader.header().tid(&args[2]).ok_or("unknown contig")?;
            let mut last = None;
            while reader.read_record(&mut raw)? {
                if raw.flags().is_unmapped() {
                    continue;
                }
                let (Some(t), Some(pos)) = (raw.tid(), raw.pos()) else {
                    continue;
                };
                let key = (t, pos);
                assert!(last.is_none_or(|previous| previous <= key), "unsorted BAM");
                last = Some(key);
                if t > tid || (t == tid && pos.as_u32() >= end) {
                    break;
                }
                if t != tid
                    || raw.flags().raw() & 0x704 != 0
                    || raw.end_pos().unwrap().as_u32() < begin
                {
                    continue;
                }
                if !raw.cigar().any(|op| op.consumes_ref() && op.len() > 0) {
                    continue;
                }
                store.push_raw(raw.as_payload(), reader.header(), &mut ())?;
            }
        }
        let retained = store.len() as i64;
        if mates {
            let stats = store.link_mates();
            println!("M\t{}\t{}", stats.pairs, stats.ambiguous_qnames);
            for (i, record) in store.records().enumerate() {
                println!("{i}\t{}", record.mate_idx().map_or(-1, |m| m.as_usize() as i64));
            }
            return Ok(());
        }
        // Seqair requires preparation (mate links and interval metadata).
        // Charge that to loading, before timing engine creation/iteration.
        let input = store.prepare_for_pileup().input;
        let load_ns = start.elapsed().as_nanos();
        let start = Instant::now();
        let mut engine = PileupEngine::new(input, span);
        let (entries, column_sum) = if dump {
            walk::<true>(&mut engine)
        } else {
            walk::<false>(&mut engine)
        };
        sum += column_sum;
        let columns_ns = start.elapsed().as_nanos();
        if dump {
            return Ok(());
        }
        let counts = black_box((records, retained, entries, sum));
        if iteration == 0 {
            expected = counts;
        } else {
            assert_eq!(counts, expected);
            elapsed[0] += read_ns;
            elapsed[1] += load_ns;
            elapsed[2] += columns_ns;
        }
        store = engine.reclaim_allocation().unwrap_or_default();
    }
    println!(
        "{} {} {} {} {} {} {}",
        elapsed[0], elapsed[1], elapsed[2], expected.0, expected.1, expected.2, expected.3
    );
    Ok(())
}

fn write_benchmark(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    use seqair::bam::BamWriterBuilder;
    use seqair::io::CompressionLevel;
    let level = CompressionLevel::try_from(args[4].parse::<i32>()?)?;
    let iterations: usize = args[5].parse()?;
    let mut reader = BamStreamReader::new(File::open(&args[2])?)?;
    let mut raw = RawBamRecord::default();
    let mut records = Vec::new();
    let mut bytes = 0;
    while reader.read_record(&mut raw)? {
        bytes += raw.as_payload().len();
        records.push(raw.clone());
    }
    let mut elapsed = 0;
    for iteration in 0..=iterations {
        let start = Instant::now();
        if args[1] == "write-sam" {
            use std::io::{BufWriter, Write};
            use seqair::sam::format::{format_sam_header, format_sam_record};
            let mut writer = BufWriter::with_capacity(65536, File::create(&args[3])?);
            let mut line = Vec::new();
            format_sam_header(reader.header(), &mut line)?;
            writer.write_all(&line)?;
            for record in &records {
                line.clear();
                format_sam_record(record, reader.header(), &mut line)?;
                line.push(b'\n');
                writer.write_all(&line)?;
            }
            writer.flush()?;
        } else {
            let mut writer = BamWriterBuilder::to_path(Path::new(&args[3]), reader.header())
                .compression_level(level).build()?;
            for record in &records { writer.write_raw(record)?; }
            writer.finish()?;
        }
        if iteration > 0 { elapsed += start.elapsed().as_nanos(); }
    }
    println!("{elapsed} {} {bytes}", records.len());
    Ok(())
}
