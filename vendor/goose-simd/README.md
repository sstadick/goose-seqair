# Shared native SIMD primitives

Unmodified `include/goosimd_expand.h` and `include/goosimd_nibbles.h` from
[goose-simd v0.3.0](https://github.com/sstadick/goose-simd/tree/ec7cbc80e5af2af8266a3facb03144bc1c83721a),
commit `ec7cbc80e5af2af8266a3facb03144bc1c83721a`. MIT license included.

These are the same implementations used by the Goose dependency, not forks.
They are vendored because standalone native C sources cannot import another
Goose package's headers by package name. Update these files and the Goose
dependency together; do not make local algorithm changes here.
