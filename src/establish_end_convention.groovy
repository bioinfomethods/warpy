import gngs.*

gngs.VCF.filter() {
    it.update { v -> {
        if (v.info.SVTYPE == "BND") {
            if (!v.info.END) {
                def m = (v.alt =~ /[\[\]]([^:]+):(\d+)[\[\]]/);
                if (m) {
                    v.info.END = m.group(2);
                } else {
                    if (v.info.END2) {
                        v.info.END = v.info.END2;
                    }
                }
            }
            if (v.info.END == "0") {
                v.info.END = "1"
            }
            if (v.info.END == 0) {
                v.info.END = 1
            }
            if (v.info.END2) {
                v.info.remove("END2")
            }
        }
    } }
}

