<script setup lang="ts">
import { ReadInfo, Segment, SegmentGroupInfo, slug } from "./segment";
import { Options } from "./options";
import { computed, ComputedRef, onMounted, ref, watchEffect } from "vue";
import * as d3 from "d3";

const props = defineProps<{
  yScale: d3.ScaleLinear<number, number, never>;
  group: SegmentGroupInfo;
  readInfo: Map<string, ReadInfo>;
  readMin: number;
  readMax: number;
  options: Options;
  colours: Map<string, string>;
}>();

const emit = defineEmits<{
  selectedSegment: [seg: Segment]
}>();



const xScale: ComputedRef<d3.ScaleLinear<number, number, never>> = computed(() => {
  return d3.scaleLinear().domain([props.group.grpMin, props.group.grpMax]).range([props.group.winMin, props.group.winMax]);
});

function y1(seg: Segment): number {
  let y = 0;
  if (seg.strand == "+") {
    y = seg.offset;
  } else {
    y = seg.offset + seg.qlen;
  }
  const info = props.readInfo.get(seg.readid);
  if (info && info.flip) {
    y = info.length - y;
  }
  return props.yScale(y);
}

function y2(seg: Segment): number {
  let y = 0;
  if (seg.strand == "+") {
    y = seg.offset + seg.qlen;
  } else {
    y = seg.offset;
  }
  const info = props.readInfo.get(seg.readid);
  if (info && info.flip) {
    y = info.length - y;
  }
  return props.yScale(y);
}

function x1(seg: Segment): number {
  return xScale.value(seg.pos);
}

function x2(seg: Segment) {
  return xScale.value(seg.pos + seg.rlen);
}

function findColour(seg: Segment) {
  return props.colours.get(seg.readid) || "black";
}

function startMark(seg: Segment) {
  const g = props.readInfo.get(seg.readid);
  if (g && seg.offset == g.begin) {
    return "url(#dotClosed)";
  } else {
    return "";
  }
}

function endMark(seg: Segment): string {
  const g = props.readInfo.get(seg.readid);
  if (g && seg.offset + seg.qlen == g.length) {
    return "url(#dotOpen)";
  } else {
    return "";
  }
}

function dashes(_seg: Segment): string {
  return "0";
}

const lineSegments = ref(null);

onMounted(() => {
  watchEffect(() => {
    d3.select(lineSegments.value)
      .selectAll("line")
      .data(props.group.segments, (seg) => slug(seg as Segment))
      .enter()
      .append("line")
      .attr("x1", (d) => x1(d))
      .attr("x2", (d) => x2(d))
      .attr("y1", (d) => y1(d))
      .attr("y2", (d) => y2(d))
      .attr("stroke", (d) => findColour(d))
      .attr("stroke-width", 2)
      .attr("marker-start", (d) => startMark(d))
      .attr("marker-end", (d) => endMark(d))
      .style("stroke-dasharray", (d) => dashes(d))
      .on("click", function (_e, seg) {
        emit('selectedSegment', seg);
      });
  });
});
</script>

<template>
  <g ref="lineSegments"> </g>
</template>
