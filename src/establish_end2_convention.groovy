import gngs.*

gngs.VCF.filter() {
    def hasEnd2 = it.header.headerLines.any { line ->
        line.startsWith("##INFO=<ID=END2,")
    };
    if (!hasEnd2) {
        it.header.addInfoHeaders(['##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of breakpoint on CHR2">']);
    }
    it.update { v -> {
        if (v.info.SVTYPE == "BND") {
            if (!v.info.END2) {
                def m = (v.alt =~ /[\[\]]([^:]+):(\d+)[\[\]]/);
                if (m) {
                    v.info.END2 = m.group(2);
                } else {
                    if (v.info.END) {
                        v.info.END2 = v.info.END;
                    }
                }
            }
            if (v.info.END2 == "0") {
                v.info.END2 = "1"
            }
            if (v.info.END2 == 0) {
                v.info.END2 = 1
            }
            if (v.info.END) {
                v.info.remove("END")
            }
        }
    } }
}
