<script setup lang="ts">
import { computedAsync } from "@vueuse/core";
//import BamAlignmentPlot from "./components/BamAlignmentPlot.vue";
import AlignmentPlot from "./components/AlignmentPlot.vue";
import { Options } from "./components/options";
import { Segment } from "./components/segment";

const options: Partial<Options> = {};

/*
const manifest = computedAsync(async () => {
  const res: Map<string, string> = new Map();
  const resp = await fetch("/manifest.json");
  const stuff = await resp.json();
  if (typeof stuff == "object") {
    const keys = Object.keys(stuff);
    for (const key of keys) {
      const value = stuff[key];
      if (typeof value == "string") {
        res.set(key, value);
      } else {
        console.log(`warning: dropping ${key} from manifest.json`);
      }
    }
  }
  return res;
});

const samples = computed(() => {
  if (manifest.value) {
    const res: string[] = [...manifest.value.keys()];
    return res;
  } else {
    const res: string[] = [];
    return res;
  }
});

const selectedSample = defineModel<string>();

const selectedBam = computed(() => {
  if (selectedSample.value && manifest.value) {
    return manifest.value.get(selectedSample.value);
  } else {
    return undefined;
  }
});

function validLocus(txt: string): string | boolean {
  const m = (txt || "").match(/^(chr([0-9]+|X|Y|MY)):([0-9]+)[-]([0-9]+)$/);
  if (!m) {
    return false;
  }
  //const chrom = m[1];
  const begin = m[3];
  const end = m[4];
  if (end != undefined) {
    if (parseInt(begin) > parseInt(end)) return "end must be greater than start";
  }
  return true;
}

function validLoci(txt: string): string | boolean {
  const parts = (txt || "").split(/[ ]+/);
  for (const part of parts) {
    const r = validLocus(part.trim());
    if (r != true) {
      return r;
    }
  }
  return true;
}

const accessToken =
  "eyJhbGciOiJSUzI1NiIsInR5cCIgOiAiSldUIiwia2lkIiA6ICIxZzNfejhDY2k5VGRXcDhNSVFMZ21OQl95NG9vZGlIWTJ5OVMwYk5OSUlvIn0.eyJleHAiOjE3MjEwMjU3NDgsImlhdCI6MTcyMTAyMzk0OCwiYXV0aF90aW1lIjoxNzE2Nzg5OTQ3LCJqdGkiOiI5NWRhODE5Ni02NmM4LTRiYzItOGUxYy1lZTdjZDU2MGZkYWEiLCJpc3MiOiJodHRwczovL2tleWNsb2FrLm1jcmkuZWR1LmF1Ojg4ODgvcmVhbG1zL2Jpb2luZm9tZXRob2RzIiwiYXVkIjpbImFyY2hpZS1uYXRpdmUiLCJhcmNoaWUiXSwic3ViIjoiNjRhMGZlNTItMjVlYi00ZmU1LWFhMmItNDVhNTk4YjY1Njc0IiwidHlwIjoiSUQiLCJhenAiOiJhcmNoaWUtbmF0aXZlIiwibm9uY2UiOiIiLCJzZXNzaW9uX3N0YXRlIjoiZTcxMGJjN2EtZTVhYi00MWJmLTllN2ItNWIxNDk2ZGI3N2ExIiwiYXRfaGFzaCI6Ii1yVjdJcGxXSm1zY19RWDJzLWhtZEEiLCJhY3IiOiIxIiwic2lkIjoiZTcxMGJjN2EtZTVhYi00MWJmLTllN2ItNWIxNDk2ZGI3N2ExIiwiYWRfZ3JvdXBzIjpbIkNFTEwxLU5HUy5ETCIsIlQtZ2l0Q0ktMDJ2LkRMIiwiTVVTQzItUk5Bc2VxLkRMIiwiTkVVUjQuREwiLCJNVVNDMi1hcmNoaWVfbWFpbi5ETCIsIm5ldXIyLWFyY2hpZV9tYWluLkRMIiwibmV1cjMtbmdzLkRMIiwiVC1naXRDSS0wMXYuREwiLCJNVVNDMi1STkFzZXEtTXVzY3VsYXJEeXN0cm9waHkuREwiLCJuZXVyNC1nZW5vbWljcy1jYXMuREwiLCJtb2xlMS1hcmNoaWVfbWFpbi5ETCIsIlBORVUxLVJTVkxhYi5ETCIsInVkcC5ETCIsIkNSQU4xLU5HUy5ETCIsIkNZVE8xLWFyY2hpZV9tYWluLkRMIiwiVC1naXRDSS0wMnYtU3Vkb2Vycy5ETCIsIkNBTkMxLWdlbm9taWMuREwiLCJULWdpdENJLTAxdi1TdWRvZXJzLkRMIiwiUE5FVTEuREwiLCJCSU9JMS1OR1MuREwiLCJuZXVyMy1hcmNoaWVfbWFpbi5ETCIsImJpb2luZi1vcHMuREwiLCJCSU9JMS5ETCIsIkNBTkMyLVJOQXNlcS5ETCIsIk1VU0MyLU5HUy5ETCIsIk5HU0RhdGFBbmFseXNlcnMuREwiLCJIRUFSMi1hcmNoaWVfbWFpbi5ETCIsIlJOQVNlcURhdGFBbmFseXNlcnMuREwiLCJDQVJEMi1hcmNoaWVfbWFpbi5ETCIsIkJMT08xLWFyY2hpZV9tYWluLkRMIiwiTkVVUjMuREwiLCJUUkFOMi0xOTA2X21vc2FpY19QTC5ETCIsIlByb2plY3RNb25pdG9yaW5nLUE4NDktMjAxNy5ETCIsIm5ldXI0LWdlbm9taWNzLkRMIiwiSFBDLUdlbm9tZXMuREwiLCJDT1JEMS1hcmNoaWVfbWFpbi5ETCIsImJpb2luZi1vcHMtdmNnc19kc3Bfc3RhZ2luZy5ETCIsIkJMT08xLU5HUy5ETCIsIlBORVUxLU5HUy5ETCIsImJtbC5ETCIsIkJJT0kxLU1HSEEuREwiLCJIUENVc2VyLkRMIiwidWRwLWdlbm9taWNfZGF0YS5ETCJdLCJyZXNvdXJjZV9hY2Nlc3MiOnsiYWNjb3VudCI6eyJyb2xlcyI6WyJtYW5hZ2UtYWNjb3VudCIsIm1hbmFnZS1hY2NvdW50LWxpbmtzIiwidmlldy1wcm9maWxlIl19fSwiZW1haWxfdmVyaWZpZWQiOnRydWUsInJlYWxtX2FjY2VzcyI6eyJyb2xlcyI6WyJvZmZsaW5lX2FjY2VzcyIsImRlZmF1bHQtcm9sZXMtYmlvaW5mb21ldGhvZHMiLCJ1bWFfYXV0aG9yaXphdGlvbiJdfSwibmFtZSI6IlRvbSBDb253YXkiLCJncm91cHMiOlsiYmlvaW5mb21ldGhvZHNfZ3JvdXAiXSwicHJlZmVycmVkX3VzZXJuYW1lIjoidG9tLmNvbndheUBtY3JpLmVkdS5hdSIsImdpdmVuX25hbWUiOiJUb20iLCJmYW1pbHlfbmFtZSI6IkNvbndheSIsImVtYWlsIjoidG9tLmNvbndheUBtY3JpLmVkdS5hdSJ9.AkSrnFlBPmguk67xdjcpJWhqWZ5y8cqp2UxPAE4Igi7SnDSLfGz4BMdxU67B-Y-jh6m1l9-k4iMmNYy_0xZYIeOs-Ls9mhXH-TH65ywQk3Er_jImZiICk5sow5STegqRsJHKoLQIJGY_pKN317NCIzcAsnGCmEMlY3kGFROnhydQILh_JwHLqkxRzoGcwYtp6KgcqhvLnn2qwYIt76mBxqrLso3aFRCY7ioRX3-0ZnXANRCt_54bYaXiHSCCtGazuFoF1VA7dYqRAMsGsGQej5_CpabBpB7fV-HBxdKSRe1UM4FYlEAV00WyZq6jbHMEeXqKWEmYj33JqhtCctVcfA";

const bamHeaders = {
  Authorisation: `Bearer ${accessToken}`,
};
*/

async function fetchFileData(file: File): Promise<{ [locus: string]: Segment[] }> {
  const blob = await file.text();
  const data = JSON.parse(blob);
  return data;
}

const selectedFile = defineModel<File | null>("selectedFile");

const rawData = computedAsync(() => {
  const file = selectedFile.value;
  if (file != undefined) {
    return fetchFileData(file);
  }
});

const currentLocus = defineModel<string>("currentLocus");
/*
const rawLoci = defineModel<string>("rawLoci", { default: "chr21:44416000-44423001" });

const loci = computed(() => {
  const res: [string, number, number][] = [];
  (rawLoci.value || "").split(/[ ]+/).forEach((s) => {
    const m = s.match(/^(chr([0-9]+|X|Y|MY)):([0-9]+)[-]([0-9]+)$/);
    if (m) {
      const chrom = m[1];
      const begin = parseInt(m[3]);
      const end = parseInt(m[4]);
      res.push([chrom, begin, end]);
    }
  });
  return res;
});
*/
</script>

<template>
  <div>
    <v-file-input v-model="selectedFile" accept="text/json"></v-file-input>
    <v-card>
      <v-tabs v-model="currentLocus">
        <v-tab v-for="(_segments, locus) in rawData" :key="locus" :value="locus">{{ locus }}</v-tab>
      </v-tabs>
      <v-tabs-window v-model="currentLocus">
        <v-tabs-window-item v-for="(segments, locus) in rawData" :key="locus" :value="locus">
          <AlignmentPlot :segments="segments" :options="options" />
        </v-tabs-window-item>
      </v-tabs-window>
    </v-card>
    <!--
      <v-select label="Sample" :items="samples" v-model="selectedSample"></v-select>
      <v-text-field
        clearable
        label="Locus/Loci"
        hint="space separated for multiple"
        placeholder="chr7:6007774-6098582"
        :rules="[validLoci]"
        v-model="rawLoci"
      ></v-text-field>
      <BamAlignmentPlot v-if="selectedBam" :bam-url="selectedBam" :bam-headers="bamHeaders" :loci="loci" :options="options" />
    -->
    <v-card> </v-card>
  </div>
</template>
