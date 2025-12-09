import gngs.*

gngs.VCF.filter() {
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
            if (v.info.END) {
                v.info.remove("END")
            }
        }
    } }
}
