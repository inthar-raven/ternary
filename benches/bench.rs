use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;
use ternary::comb::necklaces_fixed_content;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("sig 5L2m4s", |b| b.iter(|| necklaces_fixed_content(black_box(&[5, 2, 4]))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);