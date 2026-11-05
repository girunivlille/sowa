use flate2::write::GzEncoder;
use flate2::Compression as GzLevel;
use xz2::write::XzEncoder;
use zstd::stream::write::Encoder as ZstdEncoder;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use crate::reader::FastX;
use crate::qualities::{qualities_to_int, quality_quantization, int_to_qualities};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Compression {
    None,
    Gzip,
    Zstd,
    Xz,
}

impl From<String> for Compression {
    fn from(comp: String) -> Self {
        match comp.as_str() {
            "gz" => Compression::Gzip,
            "xz" => Compression::Xz,
            "zstd" => Compression::Zstd,
            _ => Compression::None
        }
    }
}

type Sink = BufWriter<Box<dyn Write>>;
pub enum OutputWriter {
    Plain(Sink),
    Gzip(GzEncoder<Sink>),
    Zstd(ZstdEncoder<'static, Sink>),
    Xz(XzEncoder<Sink>),
}

impl Write for OutputWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            OutputWriter::Plain(w) => w.write(buf),
            OutputWriter::Gzip(w) => w.write(buf),
            OutputWriter::Zstd(w) => w.write(buf),
            OutputWriter::Xz(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            OutputWriter::Plain(w) => w.flush(),
            OutputWriter::Gzip(w) => w.flush(),
            OutputWriter::Zstd(w) => w.flush(),
            OutputWriter::Xz(w) => w.flush(),
        }
    }
}

impl OutputWriter {
    /// Finishes the compressed stream (writes the trailer) and flushes the file.
    /// Must be called, otherwise the output may be truncated/invalid.
    pub fn finish(self) -> io::Result<()> {
        let mut sink = match self {
            OutputWriter::Plain(w) => w,
            OutputWriter::Gzip(w) => w.finish()?,
            OutputWriter::Zstd(w) => w.finish()?,
            OutputWriter::Xz(w) => w.finish()?,
        };
        sink.flush()
    }
}
pub fn open_for_writing(path: &Path, compression: Option<Compression>,
                        compression_level: Option<u8>) -> io::Result<OutputWriter> {
    // let file = BufWriter::new(File::create(path)?);
    // Ok(if gzip {
    //     OutputWriter::Gzip(GzEncoder::new(file, Compression::default()))
    // } else {
    //     OutputWriter::Plain(file)
    // })
    let mut compression = compression.unwrap_or(Compression::None);
    let compression_level = compression_level.unwrap_or(
        match compression {
            Compression::None => 0,
            Compression::Gzip => 6,
            Compression::Zstd => 3,
            Compression::Xz => 6
        }
    );
    let raw: Box<dyn Write>;
    if path == Path::new("-") {
        raw = Box::new(io::stdout().lock());
        compression = Compression::None;
    } else {
        raw = Box::new(File::create(path)?)
    };
    let sink = BufWriter::new(raw);

    Ok(match compression {
        Compression::None => OutputWriter::Plain(sink),
        Compression::Gzip => OutputWriter::Gzip(GzEncoder::new(sink, GzLevel::new(compression_level as u32))),
        Compression::Zstd => OutputWriter::Zstd(ZstdEncoder::new(sink, compression_level as i32)?),
        Compression::Xz => OutputWriter::Xz(XzEncoder::new(sink, compression_level as u32)),
    })
}



/// Writes one record: FASTQ if `quality` is non-empty, FASTA otherwise.
/// Works with any writer, compressed or not.
pub fn write_sequence<W: Write + ?Sized>(w: &mut W, rec: &FastX, quantization: Option<u8>) -> io::Result<()> {
    if rec.quality.is_empty() {
        writeln!(w, ">{}\n{}", rec.name, rec.sequence)
    } else {
        let mut q = qualities_to_int(&rec.quality);
        let quantiz = quantization.unwrap_or(0_u8);
        if quantiz == 0 {
            writeln!(w, "@{}\n{}\n+{}\n{}", rec.name, rec.sequence, rec.separator, rec.quality)
        } else {
            quality_quantization(&mut q, quantiz);
            writeln!(w, "@{}\n{}\n+{}\n{}", rec.name, rec.sequence, rec.separator, int_to_qualities(&q))
        }
    }
}
