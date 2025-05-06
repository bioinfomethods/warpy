BEGIN {OFS="\t"}
/^#/ {print; next} 
{
    read_support = 0

    rname_tag_start = index($8, "RNAMES=")
    if (rname_tag_start > 0) {
        rnames = substr($8, rname_tag_start + 7)
        rname_ids_end = index(rnames, ";")
        if (rname_ids_end > 0) {
            rnames = substr(rnames, 1, rname_ids_end - 1)
            gsub(/[^,]/, "", rnames)
            read_support = length(rnames) + 1
        }
        else {
            exit 1
        }
    }
    
    $8 = $8 ";SUPPORT=" read_support
    print
}
