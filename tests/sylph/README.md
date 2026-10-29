# Tiny sylph database for the test profile (D99)

`test_db.syldb` (80 KB) is a sylph 0.9.0 sketch (`sylph sketch -c 1000 -g`, in
`quay.io/biocontainers/sylph:0.9.0--ha6fb395_0`) of two public NCBI RefSeq
genomes, made on 2026-10-29:

| Accession | Organism | md5 of the .fna |
|---|---|---|
| GCF_000005845.2 | *Escherichia coli* K-12 MG1655 | 92c997bcd88e983ffdb21b2712ed3736 |
| GCF_000240185.1 | *Klebsiella pneumoniae* HS11286 | 9f1efc7218b76dbe0d128337d7726ee4 |

It only tests the plumbing (the test profile would otherwise download the
3.7 GB GTDB database). Real runs use `--sylph_db` (default: GTDB r226 c1000).
