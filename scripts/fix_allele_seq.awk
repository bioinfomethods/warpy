function convert_allele_chars(allele) {
    if (allele == ".") {
        return "N"
    }

    gsub("R", "A", allele)
    gsub("Y", "C", allele)
    gsub("K", "G", allele)
    gsub("M", "A", allele)
    gsub("S", "C", allele)
    gsub("W", "A", allele)
    gsub("B", "C", allele)
    gsub("D", "A", allele)
    gsub("H", "A", allele)
    gsub("V", "A", allele)

    return allele
}

BEGIN {FS=OFS="\t"}
/^#CHROM/ { gsub(/\t[0-9]+_/, "\t", $0); print; next }
/^#/ { print; next } 
/^chr/ {
    if (!match($4, /^<.+>$/)) {
        $4 = convert_allele_chars($4);
    }
    if (!match($5, /^<.+>$/)) {
        $5 = convert_allele_chars($5);
    }
    print;
}
