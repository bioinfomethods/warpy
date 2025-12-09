import gngs.*

gngs.VCF.filter() {
    it.update { v -> {
        if(v.info.SVTYPE in ["INS", "DEL"]) {
            v.ref = v.ref[0];
            v.alt = "<" + v.info.SVTYPE + ">";
        }
        else if (v.info.SVTYPE == "TRA") {
            v.info.SVTYPE = "BND";
            if (v.alt == "<TRA>") {
                v.alt = "<BND>";
            }
        }
    } }
}
