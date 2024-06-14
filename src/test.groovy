import gngs.*
import org.yaml.snakeyaml.*

init = {
    println "Running init"
    gngs.Region region = new gngs.Region('chr21:13267741-46296396')
    println "chr: ${region.chr}, from: ${region.from}, to: ${region.to}"
}

run() {
    init
}
