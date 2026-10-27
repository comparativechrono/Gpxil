# gpxil: genome pixel images

**One pixel per base.** A gpxil file stores a genome as pictures: every
chromosome is an image in which each pixel is one base, with separate images
(planes) for the reference sequence, a person's genotype, and per-base flags
and qualities. It is an ordinary multi-page TIFF, so any image viewer shows
the genome, any image library reads it as arrays, and a genomic interval can
be read without decoding the rest.

![HG002, a 1 kb tandem repeat on chr21, drawn three ways](docs/img/render-modes.png)

- **Lossless.** The reference (with soft-masking) and both haplotypes of one
  sample, including phase, indels, haploid regions and no-calls, can be
  reconstructed exactly.
- **Codes, not colours.** Pixels hold integer codes; colour comes from a
  palette at display time, in sRGB for screens or CMYK for print, so the same
  file makes a picture, a print or an analysis input.
- **Aligned.** Every plane, and every person's file, uses the reference's
  coordinates: the same pixel is the same genomic position.
- **Compact.** HG002's whole genome on GRCh38 (3.1 billion bases, 5.9 million
  calls) is an 888 MB file: 848 MB of reference, smaller than the FASTA.gz it
  came from, and 39 MB of genotype.

![HG002 on GRCh38, every chromosome](docs/img/karyotype.png)

## Documents

| | |
|---|---|
| [SPEC.md](SPEC.md) | the format specification, version 0.1 (draft) |
| [docs/rationale.md](docs/rationale.md) | why the format looks this way, with measurements |
| [schema/gpxil-0.1.schema.json](schema/gpxil-0.1.schema.json) | JSON Schema for the file metadata |
| [examples/](examples) | the worked example of the specification as a file, with its expected contents |

## Tools

[gpxiltools](https://github.com/comparativechrono/Gpxil-tools) is the reference implementation in
Python: it converts FASTA and VCF to gpxil, reads intervals, exports
sequences and VCF, validates files, renders PNG and CMYK TIFF pictures and
writes a print-ready "genome book" PDF with one dot per base.

```sh
pip install git+https://github.com/comparativechrono/Gpxil-tools
gpxil encode GRCh38.fa.gz -v sample.vcf.gz -o sample.gpxil.tif
gpxil render sample.gpxil.tif chr21:30,000,001-30,250,000 -o region.png --width 500
gpxil book sample.gpxil.tif -o sample-genome.pdf
```

## Background

The genome pixel image (GPxI) was introduced by Guenter Albrecht-Buehler, who
drew sequences one pixel per base, A black, C dark grey, T light grey and G
white, wrapping into rows like text so that repeats become visible patterns
([PLoS ONE 4(3): e4701, 2009](https://doi.org/10.1371/journal.pone.0004701)).
gpxil continues that idea as a file format for whole genomes and for
individual people's genomes; the `gpxi` palette reproduces his grey levels.

## Status

Draft. Version 0.1 covers one sample per file on a linear reference. Planned:
cohorts, overview pyramids, base-modification planes. Comments and issues are
welcome.

A gpxil file of a person, and every faithful picture or print of it, is that
person's genome. See section 13 of the specification.
