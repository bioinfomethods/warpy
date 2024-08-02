<script setup lang="ts">
import * as d3 from "d3";
import { Segment } from "./segment";
import { onMounted, ref, watchEffect } from "vue";
import { ColumnSpecification } from "./sortable";
import Sortable from "./Sortable.vue";

const props = defineProps<{
  segments: Segment[];
}>();

interface DisplayedSegment {
  unique: string;
  id: string;
  chrom: string;
  refStart: number;
  refEnd: number;
  strand: string;
  readStart: number;
  readEnd: number;
  refMapLen: number;
  readMapLen: number;
}

const columns: ColumnSpecification<DisplayedSegment>[] = [
  { kind: "str", title: "Chrom", field: "chrom", value: (seg: DisplayedSegment) => seg.chrom },
  { kind: "num", title: "RefStart", field: "refStart", value: (seg: DisplayedSegment) => seg.refStart },
  { kind: "num", title: "RefEnd", field: "refEnd", value: (seg: DisplayedSegment) => seg.refEnd },
  { kind: "str", title: "Strand", field: "strand", value: (seg: DisplayedSegment) => seg.strand },
  { kind: "num", title: "ReadStart", field: "readStart", value: (seg: DisplayedSegment) => seg.readStart },
  { kind: "num", title: "ReadEnd", field: "readEnd", value: (seg: DisplayedSegment) => seg.readEnd },
  { kind: "num", title: "RefMapLen", field: "refMapLen", value: (seg: DisplayedSegment) => seg.refMapLen },
  { kind: "num", title: "ReadMapLen", field: "readMapLen", value: (seg: DisplayedSegment) => seg.readMapLen },
  { kind: "str", title: "Read ID", field: "id", value: (seg: DisplayedSegment) => seg.id },
];

const items = ref<DisplayedSegment[]>([]);

onMounted(() => {
  watchEffect(() => {
    items.value = d3.map(props.segments, (seg) => {
      const res: DisplayedSegment = {
        unique: seg.id,
        id: seg.readid,
        chrom: seg.chrom,
        refStart: seg.pos,
        refEnd: seg.pos + seg.rlen,
        strand: seg.strand,
        readStart: seg.offset,
        readEnd: seg.offset + seg.qlen,
        refMapLen: seg.rlen,
        readMapLen: seg.qlen,
      };
      return res;
    });
  });
});
</script>

<template>
  <Sortable :columns="columns" :key-gen="(item: DisplayedSegment) => item.unique" v-model:data="items" height="4.5in"></Sortable>
</template>

<style>
.mono {
  font-family: monospace;
}
</style>
