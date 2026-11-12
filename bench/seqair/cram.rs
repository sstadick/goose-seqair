//! Seqair's native CRAM stream, not htslib, with optional raw-payload oracle.
use seqair::cram::stream::CramStreamReader;
use std::{fs::File, io::{BufReader, Write}, path::Path, time::Instant};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let iterations: usize = args[3].parse()?;
    let mut elapsed = 0;
    let (mut n, mut bases) = (0, 0);
    let mut dump = args.get(4).map(File::create).transpose()?;
    for i in 0..=iterations {
        let start = Instant::now();
        let mut reader = CramStreamReader::new(BufReader::new(File::open(&args[1])?))?;
        let prefix = Path::new(&args[1]).file_name().unwrap().to_str().unwrap();
        let mut decoder = reader.decoder(Path::new(&args[2]))?.name_prefix(prefix);
        let mut records = Vec::new();
        n = 0; bases = 0;
        while let Some(container) = reader.next_container()? {
            records.clear();
            decoder.decode(&container, &mut records)?;
            for record in &records {
                n += 1; bases += record.seq_len();
                if i == 0 {
                    if let Some(out) = &mut dump {
                        out.write_all(&(record.as_payload().len() as u32).to_le_bytes())?;
                        out.write_all(record.as_payload())?;
                    }
                }
            }
        }
        if i > 0 { elapsed += start.elapsed().as_nanos(); }
    }
    println!("{elapsed} {n} {bases}");
    Ok(())
}
