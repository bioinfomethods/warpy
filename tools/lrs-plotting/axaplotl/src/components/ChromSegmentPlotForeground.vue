<script setup lang="ts">
import { ReadItem, Segment, SegmentGroupInfo } from "./segment";
import { Options } from "./options";
import { computed, ComputedRef, onMounted, ref, watchEffect } from "vue";
import * as d3 from "d3";

const props = defineProps<{
  yScale: d3.ScaleLinear<number, number, never>;
  group: SegmentGroupInfo;
  reads: Map<string, ReadItem>;
  readMin: number;
  readMax: number;
  options: Options;
  colours: Map<string, string>;
}>();

const emit = defineEmits<{
  selectedSegment: [seg: Segment];
}>();

interface AugmentedSegment extends Segment {
  key: string;
}

const segments = computed<AugmentedSegment[]>(() => {
  return d3.map(props.group.segments, (seg) => {
    const info = props.reads.get(seg.readid);
    const flipped = info?.flipped == true;
    return { ...seg, key: seg.id + (flipped ? "!" : ".") };
  });
});

const boundingBoxId = computed<string>(() => {
  return `box-${props.group.winMin}-${props.group.winMax}-${props.options.height}`;
});

const xScale: ComputedRef<d3.ScaleLinear<number, number, never>> = computed(() => {
  return d3.scaleLinear().domain([props.group.grpMin, props.group.grpMax]).range([props.group.winMin, props.group.winMax]);
});

function y1(seg: Segment): number {
  const info = props.reads.get(seg.readid);
  const flipped = info?.flipped == true;
  let y = seg.strand == "+" ? seg.offset : seg.offset + seg.qlen;
  if (flipped) {
    y = info.length - y;
  }
  return props.yScale(y);
}

function y2(seg: Segment): number {
  const info = props.reads.get(seg.readid);
  const flipped = info?.flipped == true;
  let y = seg.strand == "+" ? seg.offset + seg.qlen : seg.offset;
  if (flipped) {
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

function startMark(seg: Segment): string {
  const item = props.reads.get(seg.readid);
  if (item && item.start == seg.id) {
    return "url(#dotClosed)";
  }
  return "";
}

function endMark(seg: Segment): string {
  const item = props.reads.get(seg.readid);
  if (item && item.end == seg.id) {
    return "url(#dotOpen)";
  }
  return "";
}

function dashes(_seg: Segment): string {
  return "0";
}

const lineSegments = ref(null);

onMounted(() => {
  watchEffect(() => {
    d3.select(lineSegments.value)
      .selectAll("line")
      .data(segments.value, (seg) => (seg as AugmentedSegment).key)
      .join(
        (enter) =>
          enter
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
              emit("selectedSegment", seg);
            }),
        (update) => update,
        (exit) => exit.remove()
      );
  });
});
</script>

<template>
  <g>
    <clipPath :id="boundingBoxId">
      <rect :x="group.winMin" :y="0" :width="group.winMax - group.winMin" :height="options.height"></rect>
    </clipPath>
    <g ref="lineSegments" :clip-path="'url(#' + boundingBoxId + ')'"></g>
  </g>
</template>
