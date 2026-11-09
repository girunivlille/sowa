use std::collections::{HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use simd_minimizers::{seq_hash};
use seq_hash::packed_seq;
use packed_seq::{PackedSeqVec, SeqVec};
use seq_hash::{KmerHasher, NtHasher};
mod reader;
mod writer;
mod qualities;
use clap::Parser;
use std::path::PathBuf;

#[derive(Parser, Debug)]
struct Args {
    /// Reads file to be sorted (FASTA or FASTQ)
    #[arg(long)]
    reads_in: String,

    /// Sorted output file
    #[arg(long)]
    reads_out: String,

    /// k-mer size (default = 21)
    #[arg(long, default_value_t = 21)]
    k: usize,

    /// Mean number of used k-mers per read (default = 15)
    #[arg(long, default_value_t = 15)]
    mean_nb_kmers_per_read: usize,

    /// Statistics about all k-mers (output file, default = None)
    #[arg(long)]
    nb_kmers_stats: Option<PathBuf>,

    /// Statistics about the built buckets (output file, default = None)
    #[arg(long)]
    buckets_stats: Option<PathBuf>,

    /// Statistics about the k-mers used (filtered k-mers) (output file, default = None)
    #[arg(long)]
    filtered_kmers_in_reads_stats: Option<PathBuf>,
}

fn mediane(map: &HashMap<String, usize>) -> usize {
    let mut abundances: Vec<usize> = map.values().filter(|&&v| v > 5).copied().collect();
    abundances.sort_unstable();

    let mid = abundances.len() / 2;

    if abundances.len() % 2 == 0 {
        (abundances[mid - 1] + abundances[mid]) / 2
    } else {
        abundances[mid]
    }
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();

    let fasta_file = &args.reads_in;
    let reads_out = &args.reads_out;
    let k = args.k;
    let mean_nb_kmers_per_read = args.mean_nb_kmers_per_read;

    let mut out_stats_nb_kmers: Option<BufWriter<File>> = 
    match args.nb_kmers_stats.as_ref() {
        Some(filename) => Some(BufWriter::new(File::create(filename)?)),
        None => None,
    };

    let mut out_stats_filtered_kmers_in_reads: Option<BufWriter<File>> = 
    match args.filtered_kmers_in_reads_stats.as_ref() {
        Some(filename) => Some(BufWriter::new(File::create(filename)?)),
        None => None,
    };

    let mut b_stats: Option<BufWriter<File>> = 
    match args.buckets_stats.as_ref() {
        Some(filename) => Some(BufWriter::new(File::create(filename)?)),
        None => None,
    };

    eprintln!("Step 1 : Set parameters ");

    /*  First pass on reads considering p=0.1% k-mers. 
        Count k-mers abundances.
        Define reads number, k-mers min abundance value, approximate sequencing depth,
        number of k-mers to keep after filter, percentage of k-mers to keep after filter */

    let kmers_percentage = 0.1;

    // threshold on hashing value to keep 0.1% of k-mers
    let mut t: u32 = u32::MAX / ((100.0 / kmers_percentage) as u32);

    let mut kmer_count: HashMap<String, usize> = HashMap::new();
    let mut nb_reads = 0usize;

    let mut reads_reader = reader::FastxReader::new(Path::new(&fasta_file));
    let mut record = reader::FastX::default();

    let hasher = NtHasher::<false>::new(k);

    while reads_reader.read_sequence(None, &mut record).expect("Error reading file") {
        if record.sequence.len() < k {
            continue;
        }
        let sequence = record.sequence.trim().to_ascii_uppercase();


        nb_reads += 1;

        let packed_seq = PackedSeqVec::from_ascii(sequence.as_bytes());
        let hashes: Vec<_> = hasher.hash_kmers_simd(packed_seq.as_slice(), 1).collect();

        let mut position = 0usize;

        for hashed_kmer in hashes {

            if hashed_kmer < t {
                // retrieve k-mer in the sequence if hash value < t, from hash position in the vector
                let kmer = &sequence[position..position + k];
                let entry = kmer_count.entry(kmer.to_string()).or_insert(0);
                let count = *entry + 1;
                *entry = count;
            }
            position += 1;
        }
    }


    /*  Trying to find a local minimum of abundances (k-mers with abundance below this threshold are considered errors).
        Starting from k-mers that are only present once. 
        If the number of k-mers present X times is greater than the number of k-mers present X-1 times, min_abundance = X-1 */

    let mut min_abundance = 1;
    let mut nb_min_abudance = 0usize;
    let mut sum_abundance = 0usize;
    let mut sum_abundance_inf_min = 0usize;

    for (_kmer, count) in &kmer_count{
        sum_abundance += *count;
        if *count == min_abundance{
            nb_min_abudance += 1;
        }
    }

    let mut nb_inf_min = nb_min_abudance;
    sum_abundance_inf_min += min_abundance*nb_min_abudance;

    let mut current_abundance = 0usize;
    loop{
        min_abundance += 1;
        for (_kmer, count) in &kmer_count{
            if *count == min_abundance{
                current_abundance += 1;
            }
        }
        if (current_abundance as f32) < (nb_min_abudance as f32 * 0.75) {
            nb_min_abudance = current_abundance;
            nb_inf_min += current_abundance;
            sum_abundance_inf_min += min_abundance*current_abundance
        }
        else{
            min_abundance -= 1;
            break;
        }
        current_abundance = 0;
    }

    let mut nb_sup_min_ab = kmer_count.len() - nb_inf_min;
    let mean_abundance = (sum_abundance-sum_abundance_inf_min)/(nb_sup_min_ab);
    let nb_kmers = mean_nb_kmers_per_read * nb_reads / mean_abundance;

    let mut new_p = ((nb_kmers as f32)*kmers_percentage)/(nb_sup_min_ab as f32);
    if new_p > 100.0 {
        new_p = 100.0
    }

    drop(kmer_count);

    /*  End of first pass.
        Second pass : filter k-mers greater than min_abudance using the calculated percentage */

    eprintln!("Step 2 : Filter kmers");
    
    t = u32::MAX / ((100.0 / new_p) as u32);

    let mut kmer_count: HashMap<String, usize> = HashMap::new();

    nb_sup_min_ab = 0;
    let mut reader_position = Some(0);

    while reads_reader.read_sequence(reader_position, &mut record).expect("Error reading file") {
        reader_position = None;

        let sequence = record.sequence.trim().to_ascii_uppercase();

        if sequence.len() < k {
            continue;
        }

        let packed_seq = PackedSeqVec::from_ascii(sequence.as_bytes());
        let hashes: Vec<_> = hasher.hash_kmers_simd(packed_seq.as_slice(), 1).collect();

        let mut position = 0usize;

        for hashed_kmer in hashes {

            if hashed_kmer < t {
                let kmer = &sequence[position..position + k];
                let entry = kmer_count.entry(kmer.to_string()).or_insert(0);
                let count = *entry + 1;
                *entry = count;

                if count == min_abundance+1 {
                    nb_sup_min_ab += 1;
                }
            }
            position += 1;
        }
    }

    println!("Expected number of k-mers : {}\nNumber of k-mers counted with p = {} : {}\nNumber of k-mers lower than min {} (errors) : {}\nNumber of k-mers greater than min (solid kmers) : {}", nb_kmers, new_p, kmer_count.len(), min_abundance, kmer_count.len()-nb_sup_min_ab, nb_sup_min_ab);

    let mut random_kmers: HashMap<String, isize> = HashMap::new();
    let mut good_reads: HashMap<String, Vec<usize>> = HashMap::new();
    let mut num_id = 0usize;
    let mut sum_abundances=0usize;

    for (kmer, count) in &kmer_count {
        if *count > min_abundance {
            sum_abundances += *count;
            random_kmers.insert(kmer.clone(), -1);
            good_reads.entry(kmer.clone()).or_default();
            num_id += 1;

            if num_id == nb_kmers {
                break;
            }
        }
    }

    if let Some(writer) = out_stats_nb_kmers.as_mut(){
        writeln!(writer, "total_kmers_counted\tmin_abundance\tnb_inf_min_abundance\tnb_sup_min_abundance\tnb_selected_kmers\tmean_abundance\tmedian_abundance").unwrap();
        write!(writer, "{}\t{}\t{}\t{}\t", kmer_count.len(), min_abundance, kmer_count.len()-nb_sup_min_ab, nb_sup_min_ab).unwrap();
        writeln!(writer, "{}\t{}\t{}", good_reads.len(), sum_abundances/nb_sup_min_ab, mediane(&kmer_count)).unwrap();
    }

    drop(kmer_count);

    eprintln!("kmers retrieved {}", random_kmers.len());

    let mut kmers_order: Vec<Vec<String>> = Vec::new();
    let mut ext: HashMap<usize, (String, String)> = HashMap::new();
    let mut nb_comp: usize = 0;
    let mut reads_positions: Vec<u64> = Vec::new();

    eprintln!("Step 3 : Search kmers and build graph");
    /*  Search filtered k-mers in the reads. 
        The first filtered k-mer found in each read represents the id of its bucket 
        The other ones create links in the graph of k-mers (used to order buckets in the end) */

    let mut read_nb = 0usize;

    reader_position = Some(0);
    while reads_reader.read_sequence(reader_position, &mut record).expect("Error reading file") {
        reader_position = None;
        let position = record.position;
        reads_positions.push(position);

        let sequence = record.sequence.trim().to_ascii_uppercase();

        let packed_seq = PackedSeqVec::from_ascii(sequence.as_bytes());
        let hashes: Vec<_> = hasher.hash_kmers_simd(packed_seq.as_slice(), 1).collect();

        let mut pattern_nb = 0usize;
        let mut last_kmer: Option<String> = None;

        let mut position = 0usize;
        for hashed_kmer in hashes {

            if hashed_kmer < t {
                let current_kmer = &sequence[position..position + k].to_string();

                if good_reads.contains_key(current_kmer){
                    // Place the read in its bucket
                    if pattern_nb == 0 {
                        good_reads.entry(current_kmer.to_string()).or_default().push(read_nb);

                    } else {
                        // Extend the graph
                        let last = last_kmer.as_ref().unwrap();
                        if last != current_kmer{
                            let last_comp = *random_kmers.get(last).unwrap();
                            let current_comp = *random_kmers.get(current_kmer).unwrap();

                            let last_right_ext = last_comp != -1 && ext.get(&(last_comp as usize)).map(|e| e.1 == *last).unwrap_or(false);
                            let current_left_ext = current_comp != -1 && ext.get(&(current_comp as usize)).map(|e| e.0 == *current_kmer).unwrap_or(false);

                            // Merge 2 components of the graph (last k-mer is the end of a component and current k-mer is the beggining of an other)
                            if last_right_ext && current_left_ext {
                                if last_comp != current_comp {
                                    let old_component =kmers_order[current_comp as usize].clone();

                                    for kmer in old_component {
                                        random_kmers.insert(kmer,last_comp);
                                    }

                                    let left = ext[&(last_comp as usize)].0.clone();
                                    let right = ext[&(current_comp as usize)].1.clone();

                                    ext.insert(last_comp as usize, (left, right));
                                    let mut merged = kmers_order[last_comp as usize].clone();
                                    merged.extend(kmers_order[current_comp as usize].clone());

                                    kmers_order[last_comp as usize] = merged;
                                    kmers_order[current_comp as usize].clear();
                                    ext.insert(current_comp as usize, (String::new(), String::new()));
                                }
                            }

                            // Add current k-mer at the end of a component (current k-mer is not in a component and last k-mer is at the end of a component)
                            if last_right_ext && current_comp == -1 {
                                kmers_order[last_comp as usize].push(current_kmer.clone());
                                random_kmers.insert(current_kmer.clone(), last_comp);
                                let left = ext[&(last_comp as usize)].0.clone();
                                ext.insert(last_comp as usize, (left, current_kmer.clone()));
                            }

                            // Add current k-mer at the begginning of a component
                            if current_left_ext && last_comp == -1 {
                                let mut component = vec![last.clone()];
                                component.extend(kmers_order[current_comp as usize].clone());
                                kmers_order[current_comp as usize] = component;
                                random_kmers.insert(last.clone(), current_comp);
                                let right = ext[&(current_comp as usize)].1.clone();
                                ext.insert(current_comp as usize, (last.clone(), right));
                            }

                            // Create a new component with the current and last k-mers (neither of them are in a component)
                            if last_comp == -1 && current_comp == -1 {
                                kmers_order.push(vec![last.clone(), current_kmer.clone()]);

                                random_kmers.insert(last.clone(), nb_comp as isize);
                                random_kmers.insert(current_kmer.clone(), nb_comp as isize);
                                ext.insert(nb_comp, (last.clone(), current_kmer.clone()));
                                nb_comp += 1;
                            }
                        }
                    }

                    last_kmer = Some(current_kmer.to_string());
                    pattern_nb += 1;
                }
            }
            position += 1;
        }
        if let Some(writer) = out_stats_filtered_kmers_in_reads.as_mut(){
            writeln!(writer, "{}", pattern_nb).unwrap();
        }
        read_nb += 1;
    }

    eprintln!("Step 4 : Write reads");

    let mut nb_in_buckets = 0usize;
    let mut nb_in_solo_buckets = 0usize;
    let mut nb_out_buckets = 0usize;

    let mut written_reads = vec![false; reads_positions.len()];
    let compression = if reads_reader.is_compressed() {Some(writer::Compression::Gzip)} else {Some(writer::Compression::None)};
    let mut out_fa = writer::open_for_writing(Path::new(&reads_out), compression, None).expect("Cannot create output file");

    if let Some(writer) = b_stats.as_mut(){
        writeln!(writer, "bucket_id\tnb_reads\tcomp_num").unwrap();
    }

    let mut comp_num = 0usize;

    for comp in &kmers_order {
        if !comp.is_empty() {
            comp_num += 1;
        }

        for kmer in comp {
            let reads = good_reads.get(kmer).cloned().unwrap_or_default();

            if reads.len()>0{
                if let Some(writer) = b_stats.as_mut(){
                    writeln!(writer, "{}\t{}\t{}", kmer, reads.len(), comp_num).unwrap();
                }

                for read_id in reads {
                    reads_reader.read_sequence(Some(reads_positions[read_id]), &mut record).expect("Error reading file");
                    writer::write_sequence(&mut out_fa, &record, None).expect("Cannot write to output file");

                    written_reads[read_id] = true;
                    nb_in_buckets += 1;
                }
            }
        }
    }

    for (kmer, component) in &random_kmers {
        if *component == -1 {
            let reads = good_reads.get(kmer).cloned().unwrap_or_default();

            if reads.len()>0{
                if let Some(writer) = b_stats.as_mut(){
                    writeln!(writer, "{}\t{}\t0", kmer, reads.len()).unwrap();
                }

                for read_id in reads {
                    reads_reader.read_sequence(Some(reads_positions[read_id]), &mut record).expect("Error reading file");
                    writer::write_sequence(&mut out_fa, &record, None).expect("Cannot write to output file");

                    written_reads[read_id] = true;
                    nb_in_solo_buckets += 1;
                }
            }
        }
    }

    for (read_id, written) in written_reads.iter().enumerate() {
        if !written {
            reads_reader.read_sequence(Some(reads_positions[read_id]), &mut record).expect("Error reading file");
            writer::write_sequence(&mut out_fa, &record, None).expect("Cannot write to output file");

            nb_out_buckets += 1;
        }
    }
    out_fa.finish().expect("Cannot close output file");
    
    if let Some(writer) = b_stats.as_mut(){
        writeln!(writer, "outside\t{}\t0", nb_out_buckets).unwrap();
    }
    eprintln!(
        "nb kmers : {}\nnb reads : {}\nnb selected kmers : {}\nnb reads in linked buckets : {}\nnb reads in solo buckets : {}\nnb reads in 0 bucket : {}",
        nb_kmers,
        reads_positions.len(),
        good_reads.len(),
        nb_in_buckets,
        nb_in_solo_buckets,
        nb_out_buckets
    );

    Ok(())
}