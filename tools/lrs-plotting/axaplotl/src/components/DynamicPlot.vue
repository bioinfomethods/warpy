<script setup lang="ts">
/// <reference lib="es2021" />
import { computed, ref } from "vue";
import { Options } from "./options";
import { BamFile } from "@gmod/bam";
import { BlobFile } from "generic-filehandle";
import { Locus, makeSegment, parseLocus, RawSegment, Segment } from "./segment";
import AlignmentPlot from "./AlignmentPlot.vue";
import * as d3 from "d3";
import { scanSegments } from "./scanner";
import { computeSha1 } from "./utils";

defineProps<{
  options: Partial<Options>;
}>();

function checkFileNames(bamName: string, baiName: string): boolean | string {
  if (!bamName.endsWith(".bam")) {
    return "BAM file must end with .bam";
  }
  if (!baiName.endsWith(".bai")) {
    return "BAI file must end with .bai";
  }
  if (baiName.startsWith(bamName)) {
    // foo.bam, foo.bam.bai case
    return true;
  }
  const bamBaseName = bamName.slice(0, -4);
  const baiBaseName = baiName.slice(0, -4);
  if (bamBaseName != baiBaseName) {
    return `BAM and BAI file names were not corresponding`;
  }
  return true;
}

function validFileSelection(value: File[]): boolean | string {
  if (value.length != 2) {
    return "BAM and BAI file must be selected";
  }
  const file1 = value[0];
  const file2 = value[1];
  if (file1.name.endsWith(".bam")) {
    return checkFileNames(file1.name, file2.name);
  }
  if (file2.name.endsWith(".bam")) {
    return checkFileNames(file2.name, file1.name);
  }
  return "exactly 1 BAM and 1 BAI file must be selected.";
}

const bamAndBaiFile = defineModel<File | File[] | null>("bamAndBaiFile");

const inputFiles = computed<[File, File] | undefined>(() => {
  if (bamAndBaiFile.value && Array.isArray(bamAndBaiFile.value) && validFileSelection(bamAndBaiFile.value) == true) {
    const file1 = bamAndBaiFile.value[0];
    const file2 = bamAndBaiFile.value[1];
    if (file1.name.endsWith(".bam")) {
      const res: [File, File] = [file1, file2];
      return res;
    }
    if (file2.name.endsWith(".bam")) {
      const res: [File, File] = [file2, file1];
      return res;
    }
  }
});

const bam = computed<BamFile | undefined>(() => {
  if (inputFiles.value) {
    const bamHandle = new BlobFile(inputFiles.value[0]);
    const baiHandle = new BlobFile(inputFiles.value[1]);
    const res = new BamFile({ bamFilehandle: bamHandle, baiFilehandle: baiHandle });
    return res;
  }
});

const rawLocusString = defineModel<string>("rawLociString", { default: "" });

function validLocus(txt: string): string | boolean {
  const m = (txt || "").match(/^(chr([0-9]+|X|Y|MY)):([0-9,]+)[-]([0-9,]+)$/);
  if (!m) {
    return false;
  }
  //const chrom = m[1];
  const begin = m[3].replaceAll(",", "");
  const end = m[4].replaceAll(",", "");
  if (end != undefined) {
    if (parseInt(begin) > parseInt(end)) return "end must be greater than start";
  }
  return true;
}

const locus = computed<Locus | undefined>(() => {
  if (rawLocusString.value && validLocus(rawLocusString.value) == true) {
    const res = parseLocus(rawLocusString.value.trim().replaceAll(",", ""));
    if (res) {
      return res;
    }
  }
});

const loci = computed<Locus[]>(() => {
  const res: Locus[] = [];
  if (locus.value) {
    const l = locus.value;
    const width = l.end - l.start;
    if (width <= 1000) {
      res.push(l);
    } else {
      const lhs: Locus = { chrom: l.chrom, start: d3.max([1, l.start - 10]) || 1, end: l.start + 10 };
      const rhs: Locus = { chrom: l.chrom, start: d3.max([1, l.end - 10]) || 1, end: l.end + 10 };
      res.push(lhs);
      res.push(rhs);
    }
  }
  return res;
});

const locusString = computed<string | undefined>(() => {
  if (locus.value) {
    return `${locus.value.chrom}:${locus.value.start}-${locus.value.end}`;
  }
});

const segments = ref<Segment[]>();

const scanningNow = ref<boolean>(false);

async function doScan() {
  if (bam.value && loci.value) {
    scanningNow.value = true;
    console.log("scanning...");
    const segs = await scanSegments(bam.value, loci.value);
    segments.value = d3.map(segs, makeSegment);
    scanningNow.value = false;
  }
}

const jsonBlob = computed<string>(() => {
  if (locus.value && segments.value) {
    const loc = locus.value;
    const data = segments.value;
    const res: { [locus: string]: RawSegment[] | null } = {};
    res[`${loc.chrom}:${loc.start}-${loc.end}`] = data;
    return JSON.stringify(res);
  } else {
    return JSON.stringify(null);
  }
});

async function saveJsonBlob() {
  if (locusString.value) {
    const sha = await computeSha1(locusString.value);
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(jsonBlob.value);
    const dlAnchorElem = document.getElementById("downloadAnchor");
    dlAnchorElem?.setAttribute("href", dataStr);
    dlAnchorElem?.setAttribute("download", `${sha.slice(-12)}.json`);
    dlAnchorElem?.click();
  }
}
</script>

<template>
  <v-sheet>
    <v-card>
      <v-card-text>
        <v-file-input
          v-model="bamAndBaiFile"
          label="BAM & BAI"
          hint="select a BAM and corresponding BAI file with long read alignments"
          accept=".bam,.bai"
          multiple
          :rules="[validFileSelection]"
        ></v-file-input>
        <v-text-field v-model="rawLocusString" :rules="[validLocus]" label="Locus" hint="Locus to scan" clearable></v-text-field>
        <v-btn :disabled="!locus || !bam" @click="doScan()">Scan Alignments</v-btn>
        <span :style="{ visibility: scanningNow ? 'visible' : 'hidden', 'margin-left': '1rem' }">
          <v-progress-circular indeterminate></v-progress-circular>
        </span>
      </v-card-text>
      <v-card-text v-if="locusString && segments">
        Scanning the locus {{ locusString }} yields {{ segments?.length || 0 }} aligned segments.
        <v-btn icon="mdi-download" size="x-small" @click="saveJsonBlob()"></v-btn>
        <a id="downloadAnchor" style="display: none"></a>
      </v-card-text>
    </v-card>

    <v-sheet v-if="locusString && segments">
      <AlignmentPlot :locus="locusString" :segments="segments" :options="options"></AlignmentPlot>
    </v-sheet>
  </v-sheet>
</template>
