use std::fs::{File};
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path};
use rapidgzip_core::{Decoder, IndexOptions, IndexedReader};
use std::sync::Arc;

#[derive(Debug, Default)]
pub struct FastX {
    pub name: String,
    pub quality: String,   // empty for FASTA
    pub sequence: String,
    pub separator: String, // text after '+' for FASTQ (usually empty); empty for FASTA
    pub position: u64
}

impl FastX {
    pub fn clear(&mut self) {
        self.name.clear();
        self.quality.clear();
        self.sequence.clear();
        self.separator.clear();
    }
}

/// True if the file starts with the gzip magic bytes 0x1f 0x8b.
fn is_gzip(path: &Path) -> io::Result<bool> {
    let mut magic = Vec::with_capacity(2);
    File::open(path)?.take(2).read_to_end(&mut magic)?; // handles files < 2 bytes
    Ok(magic == [0x1f, 0x8b])
}
   

pub trait ReadSeek: Read + Seek {}
impl<T: Read + Seek> ReadSeek for T {}

type FastReader = BufReader<Box<dyn ReadSeek>>;

#[derive(Debug,Clone)]
pub enum FastxType {
    FASTA, FASTQ, UNKNOWN
}
pub struct FastxReader {
    reader: FastReader,
    filetype: FastxType,
    compressed: bool
}
impl FastxReader {

    pub fn new(path: &Path) -> Self {
        let reader: Box<dyn ReadSeek>;
        let mut is_compressed: bool = false;
        if is_gzip(path).expect("Cannot determine format") {
            let source = Arc::new(File::open(path).expect("Cannot open file"));
            let decoder = Decoder::builder().decoder_threads(4).build().expect("Cannot create decoder");
            let report = decoder.decode_with_index(&source, &mut io::sink(), IndexOptions::default()).expect("Cannot index");
            reader = Box::new(IndexedReader::new(source, report.index).expect("Cannot create reader"));
            is_compressed = true;
        } else {
            reader = Box::new(File::open(path).expect("Cannot open file"))
        };        
        Self {
            reader: BufReader::new(reader),
            filetype: FastxType::UNKNOWN,
            compressed: is_compressed
        }
    }

    fn read_trimmed_line(&mut self, out: &mut String) -> io::Result<bool> {
        out.clear();
        let n: u64 = self.reader.read_line(out)? as u64;
        if n == 0 {
            return Ok(false);
        }
        while out.ends_with('\n') || out.ends_with('\r') {
            out.pop();
        }
        Ok(true)
    }

    /// Reads one FASTA or FASTQ record into `rec`.
    ///
    /// - pos: seek to that (decompressed) offset first. It must be the
    ///   start of a record. default: read at the current position.
    /// - Returns Ok(true) if a record was read, Ok(false) at end of file.
    pub fn read_sequence(
        &mut self,
        pos: Option<u64>,
        rec: &mut FastX,
    ) -> io::Result<bool> {
        if let Some(p) = pos {
            self.reader.seek(SeekFrom::Start(p))?;
        }
        rec.clear();
        rec.position = self.reader.stream_position()?;
        loop {
            if !self.read_trimmed_line(&mut rec.name)? {
                return Ok(false);
            }
            if !rec.name.is_empty() {
                break;
            }
        }
        let marker = rec.name.remove(0); // '>' or '@' (ASCII, so 1 byte)

        let mut line = String::new();
        match marker {
            '>' => {
                loop {
                    let next = self.reader.fill_buf()?;
                    if next.is_empty() || next[0] == b'>' {
                        break;
                    }
                    self.read_trimmed_line(&mut line)?;
                    rec.sequence.push_str(&line);
                }
                self.filetype = FastxType::FASTA;
            }
            '@' => {
                if !self.read_trimmed_line(&mut rec.sequence)? {
                    return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated FASTQ: no sequence"));
                }
                if !self.read_trimmed_line(&mut rec.separator)? {
                    return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated FASTQ: no '+' line"));
                }
                if !rec.separator.starts_with('+') {
                    return Err(invalid("FASTQ: third line must start with '+'"));
                }
                rec.separator.remove(0);
                if !self.read_trimmed_line(&mut rec.quality)? {
                    return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "truncated FASTQ: no quality"));
                }
                if rec.quality.len() != rec.sequence.len() {
                    return Err(invalid("FASTQ: quality and sequence lengths differ"));
                }
                self.filetype = FastxType::FASTQ;
            }
            _ => return Err(invalid("record must start with '>' (FASTA) or '@' (FASTQ)")),
        }
        Ok(true)
    }

    pub fn get_filetype(&self) -> FastxType {
        self.filetype.clone()
    }

    pub fn is_compressed(&self) -> bool {
        self.compressed
    }
} 


fn invalid(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.to_string())
}

