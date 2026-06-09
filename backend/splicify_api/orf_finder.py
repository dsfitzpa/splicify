def find_orfs(seq, min_aa_length=150):
    """Find ORFs >min_aa_length amino acids with ATG start and stop codon."""
    orfs = []
    start_codons = ['ATG']
    stop_codons = ['TAA', 'TAG', 'TGA']
    min_length = min_aa_length * 3

    def reverse_complement(dna):
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(base, base) for base in reversed(dna))

    for strand in [1, -1]:
        seq_to_search = seq if strand == 1 else reverse_complement(seq)

        for frame in range(3):
            i = frame
            while i < len(seq_to_search) - 2:
                codon = seq_to_search[i:i+3]

                if codon in start_codons:
                    j = i + 3
                    while j < len(seq_to_search) - 2:
                        stop_codon = seq_to_search[j:j+3]
                        if stop_codon in stop_codons:
                            orf_length = j + 3 - i
                            if orf_length >= min_length:
                                if strand == 1:
                                    orf_start, orf_end = i, j + 3
                                else:
                                    orf_start = len(seq) - (j + 3)
                                    orf_end = len(seq) - i

                                orfs.append({
                                    'start': orf_start,
                                    'end': orf_end,
                                    'strand': strand,
                                    'length': orf_length,
                                    'aa_length': orf_length // 3
                                })
                            break
                        j += 3
                i += 1

    # Deduplicate: keep only longest non-overlapping ORFs
    orfs = _deduplicate_orfs(orfs)
    return orfs


def _deduplicate_orfs(orfs):
    """
    Remove duplicate and overlapping ORFs, keeping the longest ones.

    Strategy:
    1. Sort by length (longest first)
    2. For each ORF, check if it significantly overlaps with any kept ORF
    3. Keep it only if overlap is <30% of the smaller ORF
    """
    if not orfs:
        return []

    # Sort by length descending (keep longest)
    sorted_orfs = sorted(orfs, key=lambda x: x['length'], reverse=True)

    kept_orfs = []
    for orf in sorted_orfs:
        # Check if this ORF significantly overlaps with any kept ORF
        has_significant_overlap = False

        for kept_orf in kept_orfs:
            # Only check overlap if on same strand
            if orf['strand'] != kept_orf['strand']:
                continue

            # Calculate overlap
            overlap_start = max(orf['start'], kept_orf['start'])
            overlap_end = min(orf['end'], kept_orf['end'])
            overlap = max(0, overlap_end - overlap_start)

            # Calculate overlap percentage relative to smaller ORF
            smaller_length = min(orf['length'], kept_orf['length'])
            if smaller_length > 0:
                overlap_pct = overlap / smaller_length

                # If >30% overlap, skip this ORF
                if overlap_pct > 0.3:
                    has_significant_overlap = True
                    break

        if not has_significant_overlap:
            kept_orfs.append(orf)

    # Sort by start position for output
    kept_orfs.sort(key=lambda x: x['start'])
    return kept_orfs
