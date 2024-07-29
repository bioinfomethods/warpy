<script setup lang="ts">
/// <reference lib="es2021" />
import { computed } from "vue";
import * as d3 from "d3";
import { Options } from "./options";
import { BamFile } from "@gmod/bam";
import { BlobFile } from "generic-filehandle";
import { Locus, parseLocus } from "./segment";
import BamAlignmentPlot from "./BamAlignmentPlot.vue";

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
    console.log(res);
    return res;
  }
});

const rawLociString = defineModel<string>("rawLociString", { default: "" });

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

const loci = computed<Locus[]>(() => {
  const res: Locus[] = [];
  if (rawLociString.value && validLoci(rawLociString.value) == true) {
    const parts = rawLociString.value.split(/[ ]+/);
    for (const part of parts) {
      const r = parseLocus(part.trim().replaceAll(",", ""));
      if (r) {
        res.push(r);
      }
    }
  }
  console.log(res);
  return res;
});
</script>

<template>
  <v-sheet>
    <v-file-input
      v-model="bamAndBaiFile"
      label="BAM & BAI"
      hint="select a BAM and corresponding BAI file with long read alignments"
      accept=".bam,.bai"
      multiple
      :rules="[validFileSelection]"
    ></v-file-input>
    <v-text-field
      v-model="rawLociString"
      :rules="[validLoci]"
      label="Loci"
      hint="white-space separated loci to scan"
      clearable
    ></v-text-field>
    <v-sheet v-if="bam">
      <BamAlignmentPlot :bam="bam" :loci="loci" :options="options"></BamAlignmentPlot>
    </v-sheet>
  </v-sheet>
</template>
