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

function convert_bnd_char(allele) {
    sub(/^R\[/, "A[", allele)
    sub(/^R\]/, "A]", allele)
    sub(/\[R$/, "[A", allele)
    sub(/\]R$/, "]A", allele)
    sub(/^Y\[/, "C[", allele)
    sub(/^Y\]/, "C]", allele)
    sub(/\[Y$/, "[C", allele)
    sub(/\]Y$/, "]C", allele)
    sub(/^K\[/, "G[", allele)
    sub(/^K\]/, "G]", allele)
    sub(/\[K$/, "[G", allele)
    sub(/\]K$/, "]G", allele)
    sub(/^M\[/, "A[", allele)
    sub(/^M\]/, "A]", allele)
    sub(/\[M$/, "[A", allele)
    sub(/\]M$/, "]A", allele)
    sub(/^S\[/, "C[", allele)
    sub(/^S\]/, "C]", allele)
    sub(/\[S$/, "[C", allele)
    sub(/\]S$/, "]C", allele)
    sub(/^W\[/, "A[", allele)
    sub(/^W\]/, "A]", allele)
    sub(/\[W$/, "[A", allele)
    sub(/\]W$/, "]A", allele)
    sub(/^B\[/, "C[", allele)
    sub(/^B\]/, "C]", allele)
    sub(/\[B$/, "[C", allele)
    sub(/\]B$/, "]C", allele)
    sub(/^D\[/, "A[", allele)
    sub(/^D\]/, "A]", allele)
    sub(/\[D$/, "[A", allele)
    sub(/\]D$/, "]A", allele)
    sub(/^H\[/, "A[", allele)
    sub(/^H\]/, "A]", allele)
    sub(/\[H$/, "[A", allele)
    sub(/\]H$/, "]A", allele)
    sub(/^V\[/, "A[", allele)
    sub(/^V\]/, "A]", allele)
    sub(/\[V$/, "[A", allele)
    sub(/\]V$/, "]A", allele)

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
        if (match($8, /SVTYPE=BND;|SVTYPE=TRA;/)) {
            $5 = convert_bnd_char($5)
        }
        else {
            $5 = convert_allele_chars($5);
        }
    }
    print;
}
