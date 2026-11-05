//rolling hashes
// main_double_hash

use std::collections::{HashMap};
use std::env;
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

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 8 {
        eprintln!(
            "Usage: {} fasta_in fasta_out k mean_nb_kmers_per_read nb_kmers_stats buckets_stats filtered_kmers_in_reads_stats",
            args[0]
        );
        std::process::exit(1);
    }

    let fasta_file = &args[1];
    let fasta_out = &args[2];
    let k: usize = args[3].parse().unwrap();
    let mean_nb_kmers_per_read: usize = args[4].parse().unwrap();
    let nb_kmers_stats = &args[5];
    let buckets_stats = &args[6];
    let filtered_kmers_in_reads_stats = &args[7];

    eprintln!("Step 1 : Set parameters ");

    let kmers_percentage = 0.1;

    let mut t: u32 = u32::MAX / ((100.0 / kmers_percentage) as u32);

    let mut kmer_count: HashMap<String, usize> = HashMap::new();
    let mut nb_reads = 0usize;
    let mut mean_reads_sizes = 0usize;

    let mut reads_reader = reader::FastxReader::new(Path::new(&fasta_file));
    let mut record = reader::FastX::default();

    let hasher = NtHasher::<false>::new(k);

    /* Première passe sur les kmers. Définit :
        le nombre de reads
        leur taille moyenne
        la valeur seuil min (min local)
        la profondeur de séquençage (~moyenne des abondances > seuil min)
        le nombre de kmers à garder (nb_reads/couverture * 10? ou qqc qui dépend de la taille des reads) 
        le nombre de kmers qui seraient gardés avec ce p (nombre total kmer_count.len - nb_inf_min)
        et donc le p à utiliser */

    while reads_reader.read_sequence(None, &mut record).expect("Error reading file") {
        if record.sequence.len() < k {
            continue;
        }
        let sequence = record.sequence.trim().to_ascii_uppercase();


        nb_reads += 1;
        mean_reads_sizes += sequence.len();

        let packed_seq = PackedSeqVec::from_ascii(sequence.as_bytes());
        let hashes: Vec<_> = hasher.hash_kmers_simd(packed_seq.as_slice(), 1).collect();

        let mut position = 0usize;

        for hashed_kmer in hashes {

            if hashed_kmer < t {
                // récupérer le kmer que si le hash est < t à partir de la position du hash dans le vecteur
                let kmer = &sequence[position..position + k];
                let entry = kmer_count.entry(kmer.to_string()).or_insert(0);
                let count = *entry + 1;
                *entry = count;
            }
            position += 1;
        }
    }

    mean_reads_sizes = mean_reads_sizes/nb_reads;

    let mut min_abundance = 1;
    let mut nb_min_abudance = 0usize;
    let mut sum_abundance = 0usize;
    let mut sum_abundance_inf_min = 0usize;
    //let mut mean_abundance = 0usize;

    for (_kmer, count) in &kmer_count{
        sum_abundance += *count;
        //mean_abundance += *count;
        if *count == min_abundance{
            nb_min_abudance += 1;
        }
    }
    let mut nb_inf_min = nb_min_abudance;
    //mean_abundance = mean_abundance/kmer_count.len();
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

    /* Deuxième passe sur les kmers : équivalente à l'ancienne première passe, filtre des kmers */
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
                // récupérer le kmer que si le hash est < t à partir de la position du hash dans le vecteur
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

    println!("Nb de kmers espéré : {}\nNb de kmers comptés avec p = {} : {}\nNb de kmers inférieurs au min {} (erreurs) : {}\nNb de kmers supérieurs au min (solid kmers) : {}", nb_kmers, new_p, kmer_count.len(), min_abundance, kmer_count.len()-nb_sup_min_ab, nb_sup_min_ab);

    let mut out_stats_nb_kmers = BufWriter::new(File::create(nb_kmers_stats).unwrap());
    writeln!(out_stats_nb_kmers, "total_kmers_counted\tmin_abundance\tnb_inf_min_abundance\tnb_sup_min_abundance\tnb_selected_kmers\tmean_abundance\tmedian_abundance").unwrap();
    write!(out_stats_nb_kmers, "{}\t{}\t{}\t{}\t", kmer_count.len(), min_abundance, kmer_count.len()-nb_sup_min_ab, nb_sup_min_ab).unwrap();

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
    writeln!(out_stats_nb_kmers, "{}\t{}\t{}", good_reads.len(), sum_abundances/nb_sup_min_ab, mediane(&kmer_count)).unwrap();

    drop(kmer_count);

    eprintln!("kmers retrieved {}", random_kmers.len());

    let mut kmers_order: Vec<Vec<String>> = Vec::new();
    let mut ext: HashMap<usize, (String, String)> = HashMap::new();
    let mut nb_comp: usize = 0;
    let mut reads_positions: Vec<u64> = Vec::new();

    eprintln!("Step 3 : Search kmers and build graph");

    let mut out_stats_filtered_kmers_in_reads = BufWriter::new(File::create(filtered_kmers_in_reads_stats).unwrap());
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

        //for mat in automaton.find_iter(&sequence) {
        let mut position = 0usize;
        for hashed_kmer in hashes {

            if hashed_kmer < t {
                let current_kmer = &sequence[position..position + k].to_string();

                if good_reads.contains_key(current_kmer){
                    if pattern_nb == 0 {
                        good_reads.entry(current_kmer.to_string()).or_default().push(read_nb);

                    } else {
                        let last = last_kmer.as_ref().unwrap();
                        if last != current_kmer{
                            let last_comp = *random_kmers.get(last).unwrap();
                            let current_comp = *random_kmers.get(current_kmer).unwrap();

                            let last_right_ext = last_comp != -1 && ext.get(&(last_comp as usize)).map(|e| e.1 == *last).unwrap_or(false);
                            let current_left_ext = current_comp != -1 && ext.get(&(current_comp as usize)).map(|e| e.0 == *current_kmer).unwrap_or(false);

                            // Fusion de deux composantes
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

                            // Ajout du kmer courant à droite
                            if last_right_ext && current_comp == -1 {
                                kmers_order[last_comp as usize].push(current_kmer.clone());
                                random_kmers.insert(current_kmer.clone(), last_comp);
                                let left = ext[&(last_comp as usize)].0.clone();
                                ext.insert(last_comp as usize, (left, current_kmer.clone()));
                            }

                            // Ajout du kmer précédent à gauche
                            if current_left_ext && last_comp == -1 {
                                let mut component = vec![last.clone()];
                                component.extend(kmers_order[current_comp as usize].clone());
                                kmers_order[current_comp as usize] = component;
                                random_kmers.insert(last.clone(), current_comp);
                                let right = ext[&(current_comp as usize)].1.clone();
                                ext.insert(current_comp as usize, (last.clone(), right));
                            }

                            // Création d'une nouvelle composante
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
        writeln!(out_stats_filtered_kmers_in_reads, "{}", pattern_nb).unwrap();
        read_nb += 1;
    }

    eprintln!("Step 4 : Write reads");

    let mut nb_in_buckets = 0usize;
    let mut nb_in_solo_buckets = 0usize;
    let mut nb_out_buckets = 0usize;

    let mut written_reads = vec![false; reads_positions.len()];
    let compression = if reads_reader.is_compressed() {Some(writer::Compression::Gzip)} else {Some(writer::Compression::None)};
    let mut out_fa = writer::open_for_writing(Path::new(&fasta_out), compression, None).expect("Cannot create output file");
    let mut b_stats =BufWriter::new(File::create(buckets_stats).unwrap());

    writeln!(b_stats, "bucket_id\tnb_reads\tcomp_num").unwrap();

    let mut comp_num = 0usize;

    for comp in &kmers_order {
        if !comp.is_empty() {
            comp_num += 1;
        }

        for kmer in comp {
            let reads = good_reads.get(kmer).cloned().unwrap_or_default();

            if reads.len()>0{
                writeln!(b_stats, "{}\t{}\t{}", kmer, reads.len(), comp_num).unwrap();

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
                writeln!(b_stats, "{}\t{}\t0", kmer, reads.len()).unwrap();

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
    
    writeln!(b_stats, "outside\t{}\t0", nb_out_buckets).unwrap();
    eprintln!(
        "nb kmers : {}\nnb reads : {}\nnb selected kmers : {}\nnb reads in linked buckets : {}\nnb reads in solo buckets : {}\nnb reads in 0 bucket : {}",
        nb_kmers,
        reads_positions.len(),
        good_reads.len(),
        nb_in_buckets,
        nb_in_solo_buckets,
        nb_out_buckets
    );
}