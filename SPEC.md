# gpxil: Genome Pixel Image format

**Version 0.1 (draft), October 2026**

gpxil stores a genome as pictures: one pixel per base, one image per
chromosome and per kind of information, in an ordinary multi-page TIFF file.
A gpxil file opens in any TIFF viewer as a picture of the genome, can be
printed, and can be read by image-analysis software as plain arrays, while
remaining a lossless representation of the reference sequence and of one
person's genotypes.

This document is the format specification. [docs/rationale.md](docs/rationale.md)
explains the design choices and the measurements behind them;
[gpxiltools](https://github.com/comparativechrono/Gpxil-tools) is the
reference implementation.

---

## Contents

1. [Conventions](#1-conventions)
2. [Overview](#2-overview)
3. [Data model](#3-data-model)
4. [Layout](#4-layout)
5. [Planes](#5-planes)
6. [Encoding VCF calls](#6-encoding-vcf-calls)
7. [Sidecar](#7-sidecar)
8. [GPXIL-TIFF container](#8-gpxil-tiff-container)
9. [Metadata](#9-metadata)
10. [Conformance](#10-conformance)
11. [Interpretation (informative)](#11-interpretation-informative)
12. [Rendering and printing (informative)](#12-rendering-and-printing-informative)
13. [Privacy and security considerations](#13-privacy-and-security-considerations)
14. [Versioning and future work](#14-versioning-and-future-work)
15. [References](#15-references)
- [Appendix A: code tables](#appendix-a-code-tables)
- [Appendix B: worked example](#appendix-b-worked-example)

---

## 1. Conventions

The key words MUST, MUST NOT, REQUIRED, SHOULD, SHOULD NOT, RECOMMENDED,
MAY and OPTIONAL are to be interpreted as described in RFC 2119 and RFC 8174
when, and only when, they appear in all capitals.

Sections marked *informative* explain or advise; everything else is
normative.

Numbers written `0x..` are hexadecimal. Bit 0 is the least significant bit.
A *nibble* is 4 bits: the *high nibble* of a byte `b` is `b >> 4`, the *low
nibble* is `b & 0x0F`.

| Term | Meaning |
|---|---|
| contig | A reference sequence: a chromosome, mitochondrion, scaffold or other named sequence. |
| offset | A 0-based position on a contig. Offset `i` is VCF `POS = i + 1` and BED start `i`. |
| plane | One per-base array over a contig, such as the reference sequence or the genotype. |
| page | One TIFF image (IFD). Each page holds one plane of one contig. |
| haplotype | One of the two copies of each chromosome in a diploid genome, numbered 1 and 2. |
| state | What a haplotype carries at a reference base: the reference base, another base, nothing, etc. |
| sidecar | Per-contig table for information that does not fit in a pixel (inserted sequences, structural variants). |
| writer / reader | Software that creates / interprets gpxil files. |

## 2. Overview

The idea goes back to Albrecht-Buehler's genome pixel images (GPxI, 2009),
in which each base of a sequence was drawn as one pixel, A black, C dark
grey, T light grey and G white, wrapping into rows like text so that
repeats line up into visible patterns. gpxil keeps the one-pixel-per-base
picture and makes it a file format:

- **Codes, not colours.** Pixels hold integer codes for bases and genotype
  states. Colours are applied when the image is shown or printed, from a
  palette; the same file can be rendered in sRGB for a screen or in CMYK for
  a press without any colour conversion of the data.
- **Planes.** Each kind of information is a separate image of the same
  shape: the reference sequence, the genotype of each base, flags, quality.
  Every plane is aligned to the reference, so the same pixel is the same
  genomic position in every plane and in every person's file.
- **An ordinary container.** A gpxil file is a valid multi-page TIFF. Pages
  are 8-bit palette images where that makes sense, so standard viewers show a
  coloured picture of the genome, and image libraries read the planes as
  arrays.
- **Random access.** Pages are stored in strips of whole rows. Any genomic
  interval is read by decoding only the strips that contain it.
- **Lossless.** The reference sequence (including soft-masking) and both
  haplotypes of the sample (substitutions, insertions, deletions, phase,
  haploid regions, no-calls) can be reconstructed exactly.

A gpxil file holds the reference and at most one sample. It is not a
replacement for VCF, CRAM or gVCF in analysis pipelines; it keeps the
haplotype content of a call set, not its record syntax or annotations
(section 6.9).

## 3. Data model

A gpxil file consists of **metadata** (section 9) and an ordered list of
**contigs**. Each contig has:

- a **name**, unique within the file;
- a **length** `L` ≥ 1, the number of bases;
- a **layout width** `W` ≥ 1 (section 4);
- one or more **planes**, each a sequence of exactly `L` values of a fixed
  data type, one value per offset `0 .. L-1`;
- an optional **sidecar** (section 7).

The standard planes are:

| Plane | Type | Content | Section |
|---|---|---|---|
| `ref` | uint8 | reference base (IUPAC code) and soft-mask bit | 5.1 |
| `gt` | uint8 | sample genotype: state of haplotype 1 and haplotype 2 | 5.2 |
| `flags` | uint8 | per-base flags: phased, insertion follows, filtered, called | 5.3 |
| `gq` | uint8 | genotype quality | 5.4 |
| `dp` | uint8 | read depth | 5.4 |
| `ps` | uint32 | phase set | 5.5 |
| `x-…` | see 5.6 | extension planes | 5.6 |

Rules:

1. A contig SHOULD have a `ref` plane. A contig without one MUST carry the
   `md5` and/or `refget` digest of its reference sequence in its metadata
   entry (section 9.3), so that readers can obtain the sequence elsewhere.
2. If a contig has a `gt` plane it MUST also have a `flags` plane.
3. `gq`, `dp` and `ps` planes MAY be present only when a `gt` plane is.
4. All contigs of a file that carry a `gt` plane describe the same sample.
5. A plane name used on several contigs MUST have the same data type on all
   of them.

## 4. Layout

Planes are stored as 2-D raster images. With layout width `W`, offset `i`
is stored at pixel

```
row = i div W      (counting from the top, starting at 0)
col = i mod W      (counting from the left, starting at 0)
```

so bases run left to right and rows run top to bottom, like lines of text.
An image has `rows = ceil(L / W)` rows and `W` columns. Pixels at
`row * W + col ≥ L` are **padding**: writers MUST set every padding pixel of
every plane to 0 and readers MUST ignore padding.

`W` is chosen by the writer and recorded per contig (section 9.3). Files
meant to be compared pixel by pixel SHOULD use the same `W`; because the
layout is a reshape of a 1-D array, a reader can always re-wrap a plane to
any other width without loss. gpxiltools uses `W = 8192` by default.

Other arrangements, such as page-sized rows for printing (section 12.4), a
width equal to a repeat period to make tandem repeats stand out as vertical
stripes, or a Hilbert curve that keeps neighbouring bases neighbouring in two
dimensions (section 11.5), are *views* computed from the raster. Version
0.1 stores only the raster layout.

## 5. Planes

Appendix A lists every code value. Values marked *reserved* MUST NOT be
written; readers SHOULD report them and MUST NOT interpret them as anything
other than "unknown".

### 5.1 `ref`: reference sequence

One byte per base:

| Bits | Meaning |
|---|---|
| 0-3 | IUPAC nucleotide code as a set over {A, C, G, T}: A = 0x1, C = 0x2, G = 0x4, T = 0x8; ambiguity codes are unions (R = A\|G = 0x5, N = 0xF). 0x0 means *no base*. |
| 4 (0x10) | soft-masked: the base was lowercase in the source FASTA (repeat-masked). |
| 5-7 | reserved, MUST be 0. |

The 4-bit set is the base encoding used for read sequences in BAM, so the
code of every IUPAC symbol is the bitwise OR of the bases it stands for
(Appendix A.1). A soft-masked base without a base (0x10) is invalid.

When reading FASTA, writers MUST map `A C G T` and the IUPAC ambiguity
letters case-insensitively as above, setting the soft-mask bit for lowercase
letters; MUST map `U` like `T`; MUST map `-` and `.` to 0x00; and MUST either
reject any other character or map it to N (0x0F) and say so in the metadata.
Writers MUST NOT drop characters: offset `i` is always base `i + 1` of the
FASTA record.

### 5.2 `gt`: genotype

One byte per base: the **high nibble** is the state of haplotype 1, the
**low nibble** the state of haplotype 2. States are relative to the
reference base at the same offset:

| State | Name | Meaning |
|---|---|---|
| 0x0 | REF | the haplotype carries the reference base |
| 0x1 | A | the haplotype carries A (a base other than the reference) |
| 0x2 | C | … C |
| 0x4 | G | … G |
| 0x8 | T | … T |
| 0xA | ABSENT | the haplotype does not exist here (one copy only: haploid region) |
| 0xC | COMPLEX | the haplotype carries something not expressible per base; see the sidecar |
| 0xD | DEL | the reference base is deleted on this haplotype |
| 0xF | UNKNOWN | the state is not known (no call) |
| others | — | reserved |

The single-bit states use the same bits as section 5.1, so the base on a
haplotype is the state itself when the state is a base, and the reference
base when the state is REF. The letters of the special states are
mnemonics: **A**bsent, **C**omplex, **D**eleted, and 0xF ("all bits") for
unknown, as N is 0xF in `ref`.

Rules:

1. A base state MUST NOT equal the reference base; a haplotype carrying the
   reference base is REF (0x0). (Where the reference base is ambiguous or N,
   any base state is allowed.)
2. For phased calls (section 5.3, PHASED) haplotype 1 and haplotype 2 are
   the first and second alleles of the phased genotype. For unphased calls,
   haplotype 1 holds the first allele written in the genotype and haplotype
   2 the second; the order then carries no biological meaning.
3. A haploid call sets haplotype 2 to ABSENT. Polyploid calls cannot be
   represented in version 0.1.
4. Insertions are not states: they happen *between* bases and are flagged in
   the `flags` plane (INS1, INS2) on the base after which they occur, with the
   inserted sequence in the sidecar.

Examples: `0x00` homozygous reference; `0x04` heterozygous with G on
haplotype 2 (`0|1`, ALT G); `0x44` homozygous G; `0x48` heterozygous with G
and T (`1|2`); `0x0D` a base deleted on haplotype 2; `0xDD` deleted on both;
`0x1A` haploid A; `0xFF` no call.

### 5.3 `flags`

One byte per base, a bit field:

| Bit | Name | Meaning |
|---|---|---|
| 0x01 | PHASED | the call covering this base was phased (`\|`), so haplotype order is meaningful |
| 0x02 | INS1 | an insertion follows this base on haplotype 1 (sequence in the sidecar) |
| 0x04 | INS2 | an insertion follows this base on haplotype 2 |
| 0x08 | FILTERED | a call covering this base did not pass its filters (VCF FILTER not PASS or `.`) |
| 0x10 | VARIANT_SITE | the base is covered by a call for this sample (any genotype) |
| 0x20-0x80 | — | reserved, MUST be 0 |

VARIANT_SITE distinguishes "called reference" (state REF, flag set) from
"reference by default" (state REF, flag clear) in files made from a VCF that
lists only variant sites.

### 5.4 `gq` and `dp`

One byte per base: genotype quality (`gq`) or read depth (`dp`) of the call
covering the base. Values above 254 MUST be stored as 254. The value 255
means *not available* and MUST be used for bases not covered by a call with
that field.

### 5.5 `ps`

One unsigned 32-bit integer per base: the phase set of the phased call
covering the base, or 0 if none. Where the VCF phase set (`FORMAT/PS`) is an
integer from 1 to 2^32−1 it MUST be stored as is. Other phase-set
identifiers MUST be replaced by distinct values not otherwise used on the
contig, and the writer SHOULD record the mapping in the contig's
`phase_set_ids` metadata. When a contig has phased calls but no `ps` plane,
all its phased calls belong to one phase set (as for phasing from a
diploid assembly or a trio).

### 5.6 Extension planes

Writers MAY add planes named `x-` followed by letters, digits, `-` or `_`,
for example `x-methylation` or `x-coverage`. An extension plane MUST have one
of the types uint8, uint16, uint32, int8, int16, int32 or float32, and MUST be
described in the `planes` metadata (section 9.2). Readers MUST ignore
extension planes they do not understand.

## 6. Encoding VCF calls

This section fixes how one sample's calls in a VCF (version 4.2 or later)
become `gt`, `flags`, `gq`, `dp`, `ps` and sidecar data, so that different
writers produce the same planes from the same input.

### 6.1 Initial state

Before any call is applied:

- every base of both haplotypes is REF;
- if the writer is given *callable regions* (for example a BED file of
  confidently called intervals), both haplotypes are UNKNOWN outside them;
- if the writer is given *haploid regions* (for example chromosome X outside
  the pseudoautosomal regions in a male), haplotype 2 is ABSENT inside them.

These initial states are not *asserted* (section 6.4). The metadata MUST say
which callable and haploid regions were applied (section 9.4).

### 6.2 Allele decomposition

A call has a REF allele of length `m` starting at offset `p = POS − 1`; it
covers the *span* `[p, p + m)`, extended to `[p, END)` for symbolic and
reference-block records that carry `INFO/END`. Each haplotype's allele is
turned into per-base states over the span as follows:

- `0` (the REF allele): REF over the span.
- `.` (missing): UNKNOWN over the span.
- `*` (spanning deletion): nothing is written; the deletion is described by
  the record that starts it.
- a sequence allele `ALT`:
  1. Trim the longest common suffix of REF and ALT, except that when the
     lengths differ the trim MUST stop before consuming the whole shorter
     allele, so the anchor base remains.
  2. Then trim the longest common prefix of what remains; let `k` be its
     length.
  3. The remaining cores `r` (from REF) and `a` (from ALT) start at offset
     `p + k`. For `j < min(|r|, |a|)`, where `r[j] ≠ a[j]` the state at
     `p + k + j` is the base `a[j]` (UNKNOWN if `a[j]` is not A, C, G or T).
  4. If `a` is longer, the surplus `a[|r|:]` is an insertion after offset
     `p + k + |r| − 1`. If that offset is before the start of the contig the
     allele is COMPLEX over the span (section 6.6).
  5. If `r` is longer, offsets `p + k + |a| .. p + k + |r| − 1` are DEL.
  6. All other offsets of the span are REF.

  Trimming the suffix first places indels as far left as the record allows,
  the same convention as VCF left-normalisation, which keeps them consistent
  with neighbouring calls inside repeats. Any alignment that yields the same
  haplotype sequence would be equally correct; fixing one makes writers agree.

### 6.3 Symbolic alleles and reference blocks

- `<NON_REF>`, `<*>` and `<X>` stand for "any other allele". A record whose
  ALT alleles are only these (a gVCF reference block) covers `[p, END)`;
  its `0` alleles assert REF over that span. A haplotype whose allele *is*
  one of them is UNKNOWN over the span.
- `<DEL…>` with `INFO/END` (or `INFO/SVLEN`): the bases `p+1 .. END−1` are DEL
  (POS is the padding base before the deletion, by VCF convention).
- Any other symbolic or breakend allele: the haplotype is COMPLEX at offset
  `p`, and a COMPLEX sidecar entry records the allele (section 7).

### 6.4 Precedence

Calls are applied in file order. Each haplotype base remembers whether a
call has *asserted* it. When a call writes a state to a base:

1. A state other than REF and UNKNOWN always replaces the current state.
   If the base already held a different state other than REF and UNKNOWN,
   this is a *conflict*: the later call wins and the writer MUST count it
   (section 9.4).
2. REF and UNKNOWN replace the current state only if no call has asserted
   the base yet.
3. The base is now asserted.

So a record whose REF allele overlaps another record's variant never
undoes that variant, whatever the record order, and overlapping records
that use `*` combine correctly. Split multi-allelic records (two records at
one position, each `0/1`) cannot be resolved and show up as conflicts.

Insertions are added at their anchor base on their haplotype; two different
insertions at the same base on the same haplotype are a conflict (last
wins).

### 6.5 Ploidy, phase and flags

- A genotype with one allele (`1`, `0`, `.`) is haploid: haplotype 2 is
  ABSENT over the span. Genotypes with more than two alleles MUST be skipped
  and counted.
- A genotype whose separators are all `|` is phased: PHASED is set over the
  span (not for haploid calls).
- VARIANT_SITE is set over the span of every call; FILTERED over the span of
  every call whose FILTER is neither `PASS` nor `.`.
- `gq` and `dp` take `FORMAT/GQ` and `FORMAT/DP` over the span (capped at
  254); `ps` takes `FORMAT/PS` over the span of phased calls.

### 6.6 COMPLEX

A haplotype is COMPLEX where its allele cannot be expressed per base: symbolic
alleles other than deletions (section 6.3), and insertions with no anchor
base. Every COMPLEX run starts with a COMPLEX sidecar entry (section 7).

### 6.7 Reference checks

Writers MUST check each REF allele against the reference (a REF base matches
when it is contained in the reference's IUPAC set, or either is N) and MUST
count mismatches in the metadata. They MAY reject files with mismatches.

### 6.8 Contig names

Writers MAY match VCF contig names to reference names that differ only by a
`chr` prefix, or by `MT`/`M`/`chrM`, and SHOULD record the alternative names
as `aliases` (section 9.3).

### 6.9 What is and is not kept

Kept exactly: the reference sequence and soft-masking; the sequence of both
haplotypes (every substitution, deletion and insertion); phasing; haploid
calls; no-calls; the FILTER pass/fail status, genotype quality, depth and
phase set per base.

Not kept: VCF record boundaries and allele syntax, IDs, QUAL, INFO, other
FORMAT fields, and other samples. Converting a gpxil file back to VCF
(section 11.4) yields records that are equivalent in haplotype content to
the input, not a copy of it.

## 7. Sidecar

Each contig with a `gt` plane MAY have a sidecar: a table of events that do
not fit in pixels. It is UTF-8 text with tab-separated columns:

```
#gpxil-sidecar	0.1
offset	hap	kind	data
5217935	1	INS	AAAAAGA
5217935	2	INS	AAAAAGA
9812201	2	COMPLEX	REF=A;ALT=<INV>;END=9813300
```

| Column | Content |
|---|---|
| `offset` | 0-based offset of the base the entry belongs to |
| `hap` | 1 or 2 |
| `kind` | `INS` or `COMPLEX` |
| `data` | for `INS`: the inserted bases, uppercase `A C G T N`; for `COMPLEX`: `;`-separated `KEY=value` items, at least `REF` and `ALT`, and `END` (VCF coordinates) when known |

Rules:

1. The first line is `#gpxil-sidecar` and the version, the second the
   column header; then one entry per line, sorted by offset, then hap, then
   kind.
2. For every base with INS1 (INS2) set there MUST be exactly one `INS` entry
   with that offset and hap 1 (2), and every `INS` entry MUST have its flag set.
3. Readers MUST ignore entries of unknown kinds and lines starting with `#`.

## 8. GPXIL-TIFF container

### 8.1 File

A gpxil file is a TIFF file (TIFF 6.0, or BigTIFF when larger than 4 GiB) in
either byte order. Writers MUST use BigTIFF when the file may exceed 4 GiB
and SHOULD use classic TIFF otherwise, because more software reads it.

The RECOMMENDED file name extension is `.gpxil.tif`, which tells people what
the file is while letting ordinary software open it; the media type is
`image/tiff`. Readers MUST recognise a gpxil file by its metadata (section
9), not by its name.

### 8.2 Pages

Every plane of every contig is one page (IFD). Pages appear contig by contig
in the order of the `contigs` metadata list; within a contig, in the order
`ref`, `gt`, `flags`, `gq`, `dp`, `ps`, then extension planes in lexical
order. Page 0 MUST be a data page and MUST hold the file metadata
(section 8.4). Readers MUST locate pages through the metadata and MUST ignore
pages it does not list.

### 8.3 Tags

Each data page has these tags:

| Tag | Value |
|---|---|
| ImageWidth (256) | `W` |
| ImageLength (257) | `ceil(L / W)` |
| BitsPerSample (258) | 8 for uint8/int8 planes, 16 or 32 according to type |
| SamplesPerPixel (277) | 1 |
| SampleFormat (339) | 1 (unsigned), 2 (signed) or 3 (float) according to type; may be omitted for unsigned |
| PhotometricInterpretation (262) | 3 (palette) for `ref`, `gt` and `flags`; 1 (BlackIsZero) for all others |
| ColorMap (320) | for palette pages: display colours (section 8.5) |
| Compression (259) | see below |
| Predictor (317) | 1 (none) or absent for `ref`, `gt` and `flags`; any for extension planes |
| ImageDescription (270) | page metadata (section 8.4) |
| PageName (285) | RECOMMENDED: `"<contig> <plane>"` for display |
| StripOffsets, RowsPerStrip, StripByteCounts or the tile tags | data layout |
| 65400 (GPXIL_SIDECAR) | on the `gt` page only, when the contig has a sidecar: the sidecar text compressed with zlib (RFC 1950), type UNDEFINED |

Tag 65400 lies in the TIFF range reserved for reusable private tags.

**Compression.** Readers MUST support no compression (1) and Deflate
(8, and the obsolete code 32946). Readers SHOULD support LZW (5),
PackBits (32773), Zstandard (50000) and LZMA (34925). Writers SHOULD use
Deflate, which every TIFF reader supports; Zstandard and LZMA give files
5-10 % smaller at the cost of compatibility with some viewers.

**Strips.** Writers SHOULD store pages in strips (whole rows), sized so that
a strip holds between 64 KiB and 1 MiB of uncompressed data; gpxiltools uses
256 KiB. Then a genomic interval maps to a run of rows and only the strips
holding those rows need decoding. Tiles are allowed; readers MUST support
both.

### 8.4 Page metadata

The ImageDescription of every page except page 0 MUST be a JSON object with
at least:

```json
{"gpxil": "0.1", "contig": "chr21", "plane": "gt", "length": 46709983}
```

The ImageDescription of page 0 MUST be the file metadata (section 9), which
also contains `contig`, `plane` and `length` for page 0 itself. JSON MUST be
written as 7-bit ASCII, with any other characters escaped (`\uXXXX`), so
that it is also a valid TIFF ASCII field.

### 8.5 Colour maps

The ColorMap of a palette page is a viewing aid only: it makes the page a
recognisable picture in an ordinary image viewer. Readers MUST NOT use it to
interpret pixel values. The `ref` colour map SHOULD follow one of the
palettes in section 12.2 and its name SHOULD be given in the `palette`
metadata key. The colours of the `gt` and `flags` maps are not specified.

## 9. Metadata

The file metadata is a JSON object (RFC 8259). [schema/gpxil-0.1.schema.json](schema/gpxil-0.1.schema.json)
is a JSON Schema for it. Readers MUST ignore keys they do not know.

### 9.1 Top level

| Key | Required | Content |
|---|---|---|
| `gpxil` | yes | format version, `"0.1"` |
| `layout` | yes | `{"type": "raster", "width": W}`; `width` is the default for contigs |
| `contigs` | yes | list of contig entries (9.3), in page order |
| `planes` | yes if extension planes are present | plane descriptions (9.2) |
| `reference` | RECOMMENDED | `{"name": "GRCh38", "source": "<file name>", ...}` |
| `sample` | yes if any `gt` plane is present, else `null` | sample description (9.4) |
| `palette` | no | name of the palette used for the `ref` colour maps |
| `generator`, `created` | no | software and ISO 8601 creation time |
| `contig`, `plane`, `length` | yes | page 0's own page metadata (8.4) |

### 9.2 Planes

`planes` maps plane names to objects with `dtype` (`"uint8"`, `"uint16"`,
`"uint32"`, `"int8"`, `"int16"`, `"int32"` or `"float32"`) and optionally
`missing` (the value meaning not available), `description` and `units`.
Extension planes MUST be described here.

### 9.3 Contig entries

| Key | Required | Content |
|---|---|---|
| `name` | yes | contig name: non-empty, no whitespace or control characters; SHOULD follow the SAM naming rules |
| `length` | yes | `L` |
| `pages` | yes | map from plane name to page index |
| `width` | no | layout width `W` of this contig, if different from `layout.width` |
| `md5` | RECOMMENDED | MD5 of the uppercase sequence, as in SAM `@SQ M5` |
| `refget` | RECOMMENDED | GA4GH refget identifier `SQ.<sha512t24u>` of the uppercase sequence |
| `aliases` | no | other names of this contig (e.g. `"21"` for `"chr21"`) |
| `ploidy` | no | `2`, `1` or `"mixed"`, as encoded |
| `vcf_contig` | no | the VCF contig name that was encoded |
| `encoding` | no | conversion counts (9.4) |
| `callable_bases` | no | number of bases inside callable regions |
| `sidecar_entries` | no | number of sidecar entries |
| `phase_set_ids` | no | map from `ps` values to original phase-set identifiers (5.5) |

### 9.4 Sample and conversion record

`sample` is `null` for a reference-only file, otherwise an object with
`id` (REQUIRED, the VCF sample name) and, as applicable, `source` (input
file name), `callable_regions` and `haploid_regions` (names of the region
files applied in section 6.1, or `null`) and `sex`.

Each contig's `encoding` object SHOULD report `records`, `conflicts`
(section 6.4), `ref_mismatches` (6.7), `polyploid_skipped` (6.5),
`skipped_out_of_range`, `spanning_deletion_alleles` and `symbolic_alleles`.

## 10. Conformance

A **conforming writer** produces files in which:

1. the TIFF structure, page order, tags and page metadata follow section 8;
2. every plane has the dimensions of section 4 and zero padding;
3. every code is valid (section 5, Appendix A): `ref` bits 5-7 clear and no
   0x10; `gt` states from the table; reserved `flags` bits clear; a base state
   never equals the reference base;
4. INS flags and sidecar `INS` entries correspond one to one (section 7);
5. the metadata contains the required keys (section 9) and, where present,
   digests that match the `ref` plane;
6. calls from a VCF are converted as in section 6.

A **conforming reader**:

1. reads the metadata from page 0 and finds pages through `contigs[].pages`;
2. decodes uncompressed and Deflate pages, in strips or tiles;
3. ignores padding, unknown metadata keys, unknown extension planes, unlisted
   pages and unknown sidecar kinds;
4. treats reserved `gt` states as UNKNOWN;
5. does not use colour maps to interpret data.

`gpxil validate` in gpxiltools checks files against this section.

## 11. Interpretation (informative)

### 11.1 Haplotype sequence

The sequence of haplotype `h` over offsets `[s, e)` is built base by base:
REF gives the reference base; A, C, G or T gives that base; UNKNOWN and
COMPLEX give N; DEL and ABSENT give nothing; after each base with INSh set,
the sidecar's inserted bases follow. Soft-mask case can be carried over from
the reference (the convention of `bcftools consensus`).

### 11.2 Genotype as one code

The union of the two haplotype bases at a position (`haplotype 1 base |
haplotype 2 base`, with DEL and ABSENT contributing nothing and UNKNOWN or
COMPLEX contributing N) is an IUPAC code: a heterozygous A/G site is R.
This is the natural single-code picture of a diploid genome, and what
rendering uses for heterozygous colours (section 12.3).

### 11.3 Dosage

Counting, at every position, how many haplotypes carry A, C, G and T gives
an `(L, 4)` array of allele dosages: the one-hot encoding used by DNA
sequence models, multiplied by ploidy. It is computed from `ref` and `gt`
and is not stored; storing it instead would lose phase, insertions and the
distinction between deleted and unknown bases (see docs/rationale.md).

### 11.4 Back to VCF

A gpxil file converts back to VCF by emitting one record per *event window*:
each substitution is a one-base window; a run of deleted bases takes the
anchor base before it; an insertion takes its anchor base; windows that
share a base are merged. REF is the reference over the window and each
haplotype contributes its sequence over the window as an allele. Runs of
UNKNOWN outside windows can be written as gVCF-style blocks
(`ALT=<*>`, `INFO/END`, `GT=./.`) so that re-encoding restores them.

### 11.5 Image analysis

Because planes are aligned to the reference, image operations are genome
operations:

- **comparison**: the pixel-wise difference of two people's `gt` planes lists
  where their genomes differ;
- **cohorts**: the mean of many dosage images is an allele-frequency image;
- **density**: a box filter over the variant mask of `gt` is a variant-density
  map, and over the G/C bits of `ref` a GC-content (isochore) map;
- **periodicity**: re-wrapping `ref` to a width equal to a repeat period turns
  tandem repeats into vertical stripes, which Albrecht-Buehler used to find
  tetra-periodic motifs;
- **two-dimensional locality**: a Hilbert-curve view keeps runs of bases
  together in 2-D blocks (a `2^k × 2^k` block is `4^k` consecutive bases), so
  convolutional filters see genomic neighbourhoods.

## 12. Rendering and printing (informative)

### 12.1 Data are not colours

A pixel value is a code; a palette turns codes into colours when the image is
displayed or printed. This separation is what lets one file serve both the
screen (sRGB) and the press (CMYK). Converting colour pixel data between
colour spaces is lossy: CMYK has four channels and RGB three, so many CMYK
values print the same colour, and ICC conversions apply gamut mapping and
black generation. Data stored as colours do not survive such conversions.
Instead, a print rendering maps codes directly to ink values chosen in CMYK,
and a screen rendering maps the same codes directly to sRGB.

### 12.2 Palettes

A palette gives the four bases, N and "no base" in both sRGB and CMYK (ink
coverage in percent); other codes are derived:

- **ambiguity codes and heterozygous sites** mix their component bases,
  either as the mean of the colours, or by *overprinting* (every component
  ink at full strength, multiplying transmittance on screen);
- **soft-masked bases** are faded towards paper white by a fixed fraction.

| Palette | A | C | G | T | Mix | Notes |
|---|---|---|---|---|---|---|
| `classic` | red | blue | green | yellow | mean | colours of the original gpxil notebooks |
| `gpxi` | black | dark grey | white | light grey | mean | Albrecht-Buehler (2009) grey levels; no fade |
| `okabe-ito` | bluish green | blue | orange | vermillion | mean | safe for colour-blind readers; IGV-like order |
| `process` | cyan | magenta | yellow | black | overprint | one process ink per base |

Exact values are in Appendix A.4. Under the `process` palette the 4-bit
IUPAC code *is* a CMYK pixel at one bit per ink: bit A lays down cyan, C
magenta, G yellow, T black. A heterozygous A/G site prints cyan and yellow;
four channels for four bases, with each code mapping to exactly one ink
combination.

### 12.3 Render modes

- **reference**: the `ref` plane.
- **genotype**: the union code of section 11.2, so homozygous bases show their
  base colour, heterozygous sites the mixture, and bases deleted on both
  haplotypes show as paper.
- **haplotypes**: two pixel rows per row of bases, haplotype 1 above
  haplotype 2. Lossless for substitutions and deletions.
- **variants**: the reference faded, with variant sites coloured by
  zygosity (heterozygous, homozygous, deletion, insertion, no call).

### 12.4 A genome book

Printed at one dot per base, an A4 page with 15 mm margins at 300 dots per
inch holds 2,126 bases per row and 6,314,220 bases per page, so the 25
sequences of GRCh38's chromosomes 1-22, X, Y and the mitochondrion (3.09
billion bases) fill 505 pages; at 600 dots per inch they fill 137. Each page
is one indexed image whose colour table is in DeviceCMYK, with interpolation
off so every base stays a square dot. Colours given in device CMYK should be
proofed with the printer's ICC profile.

## 13. Privacy and security considerations

A gpxil file with a `gt` plane is a person's genome, and so is every faithful
rendering of it. A genome book printed at one dot per base is the genome in
printed form: scanned at sufficient resolution it can be read back. Treat
files, renderings and books of individuals with the same care as the VCF or
CRAM they came from: they are identifying, and they reveal information about
relatives. Pictures intended for display or sharing should use a
reference-only or deliberately reduced rendering. Metadata may contain
sample identifiers and file names.

Readers parse JSON, a zlib stream and TIFF structures from untrusted input
and should bound memory use (a page header can claim a very large image) and
reject inconsistent sizes.

## 14. Versioning and future work

The version is `major.minor`. Minor versions add optional features (new
planes, metadata keys, sidecar kinds, state codes taken from the reserved
values) that older readers can ignore safely. A new major version is not
readable by older readers. Readers MUST reject files with a major version
they do not implement.

Planned or possible extensions:

- multi-sample (cohort) files, as stacks of `gt` planes, and a Zarr binding
  for analysis at cohort scale;
- reduced-resolution overview images (TIFF SubIFDs) for zooming viewers;
- a standard plane for base modifications (e.g. 5mC methylation from
  long-read sequencing);
- an OME-TIFF compatible profile for bioimaging software;
- polyploid genotypes and richer structural-variant support;
- preview pages for ordinary viewers.

## 15. References

- Albrecht-Buehler G (2009). Outline of a genome navigation system based on
  the properties of GA-sequences and their flanks. *PLoS ONE* 4(3): e4701.
  https://doi.org/10.1371/journal.pone.0004701
- Adobe (1992). *TIFF Revision 6.0*; and the BigTIFF extension (2007).
- SAMtools / hts-specs. *The Variant Call Format Specification*, v4.3 and v4.4.
- SAMtools / hts-specs. *Sequence Alignment/Map Format Specification*
  (BAM 4-bit base encoding).
- GA4GH. *refget* sequence retrieval specification (sha512t24u identifiers).
- Cornish-Bowden A (1985). Nomenclature for incompletely specified bases in
  nucleic acid sequences (IUPAC-IUB). *Nucleic Acids Research* 13: 3021-3030.
- Okabe M, Ito K (2002). *Color Universal Design*.
- Wagner J et al. Genome in a Bottle HG002 benchmark sets (used for the
  measurements in docs/rationale.md).

---

## Appendix A: code tables

### A.1 `ref` low nibble

| Code | Symbol | Bases | | Code | Symbol | Bases |
|---|---|---|---|---|---|---|
| 0x0 | `-` | none | | 0x8 | T | T |
| 0x1 | A | A | | 0x9 | W | A T |
| 0x2 | C | C | | 0xA | Y | C T |
| 0x3 | M | A C | | 0xB | H | A C T |
| 0x4 | G | G | | 0xC | K | G T |
| 0x5 | R | A G | | 0xD | D | A G T |
| 0x6 | S | C G | | 0xE | B | C G T |
| 0x7 | V | A C G | | 0xF | N | A C G T |

Add 0x10 for a soft-masked (lowercase) base. Valid `ref` bytes are
0x00-0x0F and 0x11-0x1F.

### A.2 `gt` nibble states

| Nibble | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | A | B | C | D | E | F |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| State | REF | A | C | – | G | – | – | – | T | – | ABSENT | – | COMPLEX | DEL | – | UNKNOWN |

(– = reserved)

### A.3 `flags` bits

`0x01` PHASED, `0x02` INS1, `0x04` INS2, `0x08` FILTERED,
`0x10` VARIANT_SITE, `0x20 0x40 0x80` reserved.

### A.4 Palette values

sRGB as 0-255 triples; CMYK as ink percentages.

| Palette | Base | sRGB | CMYK % |
|---|---|---|---|
| classic | A | 255 0 0 | 0 100 100 0 |
| | C | 0 0 255 | 100 80 0 0 |
| | G | 0 255 0 | 80 0 100 0 |
| | T | 255 255 0 | 0 0 100 0 |
| | N | 128 128 128 | 0 0 0 50 |
| gpxi | A | 0 0 0 | 0 0 0 100 |
| | C | 85 85 85 | 0 0 0 67 |
| | G | 255 255 255 | 0 0 0 0 |
| | T | 170 170 170 | 0 0 0 33 |
| | N | 170 200 235 | 30 10 0 0 |
| okabe-ito | A | 0 158 115 | 97 0 75 0 |
| | C | 0 114 178 | 100 50 0 0 |
| | G | 230 159 0 | 0 50 100 0 |
| | T | 213 94 0 | 0 80 100 0 |
| | N | 153 153 153 | 0 0 0 40 |
| process | A | 0 174 239 | 100 0 0 0 |
| | C | 236 0 140 | 0 100 0 0 |
| | G | 255 242 0 | 0 0 100 0 |
| | T | 35 31 32 | 0 0 0 100 |
| | N | 160 160 160 | 0 0 0 35 |

"No base" is paper: sRGB 255 255 255, no ink (the `gpxi` palette tints it
255 250 230 / 0 2 10 0 so that it differs from G). Soft-masked bases are faded
45 % towards paper in sRGB and their ink reduced by 45 % in CMYK (`gpxi`: no
fade). sRGB values of the `process` palette are screen approximations of the
inks.

## Appendix B: worked example

Reference `chrT` = `ACGTACGTac` (10 bases, the last two soft-masked),
layout width 4, and one sample with these calls:

```
#CHROM POS REF  ALT   FILTER GT
chrT   2   C    T     PASS   0|1
chrT   4   TAC  T     PASS   1|0
chrT   8   T    TGG   PASS   1|1
chrT   9   A    G     LowQ   ./1
```

Planes (offsets 0-9; padding offsets 10 and 11 are 0):

| offset | base | `ref` | `gt` | meaning | `flags` |
|---|---|---|---|---|---|
| 0 | A | 0x01 | 0x00 | REF\|REF | 0x00 |
| 1 | C | 0x02 | 0x08 | REF\|T | 0x11 PHASED, VARIANT_SITE |
| 2 | G | 0x04 | 0x00 | | 0x00 |
| 3 | T | 0x08 | 0x00 | anchor of the deletion | 0x11 |
| 4 | A | 0x01 | 0xD0 | DEL\|REF | 0x11 |
| 5 | C | 0x02 | 0xD0 | DEL\|REF | 0x11 |
| 6 | G | 0x04 | 0x00 | | 0x00 |
| 7 | T | 0x08 | 0x00 | insertion follows on both | 0x17 PHASED, INS1, INS2, VARIANT_SITE |
| 8 | a | 0x11 | 0xF4 | UNKNOWN\|G | 0x18 FILTERED, VARIANT_SITE |
| 9 | c | 0x12 | 0x00 | | 0x00 |

Sidecar:

```
#gpxil-sidecar	0.1
offset	hap	kind	data
7	1	INS	GG
7	2	INS	GG
```

Following section 11.1, haplotype 1 is `A C G T` (offsets 4 and 5 deleted)
`G T`, the inserted `GG`, then `N` (unknown) and `C`: **`ACGTGTGGNC`**.
Haplotype 2 is **`ATGTACGTGGGC`**. Carrying the reference's soft-mask case
over, they read `ACGTGTGGnc` and `ATGTACGTGGgc`.

As images with `W = 4`, each plane is 3 rows of 4 pixels; the `ref` page
shows `ACGT / ACGT / ac` followed by two padding pixels.

The file `examples/worked-example.gpxil.tif` in this repository contains
exactly this example.
