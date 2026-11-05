
/// Transform a quality string as in a FASTQ into a vector of integers (corresponding to the ASCII codes)
pub fn qualities_to_int(qualities: &String) -> Vec<u8> {
    let mut qual: Vec<u8> = Vec::with_capacity(qualities.len());
    for c in qualities.chars() {
        qual.push(c.to_ascii_uppercase() as u8);
    }
    return qual;
}

/// Perform the opposite compared to qualities_to_int
pub fn int_to_qualities(qualities: &Vec<u8>) -> String {
    let mut s = String::with_capacity(qualities.len());
    for q in qualities.iter() {
        s.push((*q as u32).try_into().unwrap());
    }
    return s;
}

/// Perform a quantization of the values provided in qualities
/// A quantization of x means that if n values were possible in
/// total, now n / x will be used
pub fn quality_quantization(qualities: &mut Vec<u8>, quantization_factor: u8) {
    for q in qualities.iter_mut() {
        *q = (*q as f32 / quantization_factor as f32).round() as u8 * quantization_factor;
    }
}

