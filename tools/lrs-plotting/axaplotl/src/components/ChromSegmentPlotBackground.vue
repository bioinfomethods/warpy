<script setup lang="ts">
import { flipSegmentIfNecessary, ReadItem, Segment, SegmentGroupInfo, slug } from "./segment";
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

const segments = computed<Segment[]>(() => {
  const res = d3.map(props.group.segments, (seg) => {
    const item = props.reads.get(seg.readid);
    if (item) {
      return flipSegmentIfNecessary(item, seg);
    } else {
      return seg;
    }
  });
  return res;
});

const boundingBoxId = computed<string>(() => {
  return `box-${props.group.winMin}-${props.group.winMax}-${props.options.height}`;
});

const xScale: ComputedRef<d3.ScaleLinear<number, number, never>> = computed(() => {
  return d3.scaleLinear().domain([props.group.grpMin, props.group.grpMax]).range([props.group.winMin, props.group.winMax]);
});

const panelHeight = computed(() => {
  return props.options.height - (props.options.marginBottom + props.options.marginTop);
});

const axisTransform = computed(() => {
  return `translate(0,${props.options.height - props.options.marginBottom})`;
});

function y1(seg: Segment): number {
  const y = seg.offset + (seg.strand == "+" ? 0 : seg.qlen);
  return props.yScale(y);
}

function y2(seg: Segment): number {
  const y = seg.offset + (seg.strand == "+" ? seg.qlen : 0);
  return props.yScale(y);
}

function x1(seg: Segment): number {
  return xScale.value(seg.pos);
}

function x2(seg: Segment) {
  return xScale.value(seg.pos + seg.rlen);
}

const xAxis = ref(null);
const guidesV1 = ref(null);
const guidesV2 = ref(null);
const guidesH1 = ref(null);
const guidesH2 = ref(null);

onMounted(() => {
  watchEffect(() => {
    d3.select(xAxis.value)
      .call(d3.axisBottom(xScale.value).ticks(props.group.ticks) as any)
      .selectAll("text")
      .attr("y", 0)
      .attr("x", 9)
      .attr("dy", ".35em")
      .attr("transform", "rotate(45)")
      .style("text-anchor", "start");

    d3.select(guidesV1.value)
      .selectAll("line")
      .data(segments.value, (seg) => slug(seg as Segment))
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("x1", (d) => x1(d))
            .attr("x2", (d) => x1(d))
            .attr("y1", (_d) => props.yScale(props.readMin))
            .attr("y2", (_d) => props.yScale(props.readMax))
            .attr("stroke", (_d) => props.options.guideColour)
            .attr("stroke-width", 1),
        (update) => update,
        (exit) => exit.remove()
      );
    d3.select(guidesV2.value)
      .selectAll("line")
      .data(segments.value, (seg) => slug(seg as Segment))
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("x1", (d) => x2(d))
            .attr("x2", (d) => x2(d))
            .attr("y1", (_d) => props.yScale(props.readMin))
            .attr("y2", (_d) => props.yScale(props.readMax))
            .attr("stroke", (_d) => props.options.guideColour)
            .attr("stroke-width", 1),
        (update) => update,
        (exit) => exit.remove()
      );

    d3.select(guidesH1.value)
      .selectAll("line")
      .data(segments.value, (seg) => slug(seg as Segment))
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("x1", (_d) => props.options.marginLeft)
            .attr("x2", (_d) => props.options.width - props.options.marginRight)
            .attr("y1", (d) => y1(d))
            .attr("y2", (d) => y1(d))
            .attr("stroke", (_d) => props.options.guideColour)
            .attr("stroke-width", 1),
        (update) => update,
        (exit) => exit.remove()
      );
    d3.select(guidesH2.value)
      .selectAll("line")
      .data(segments.value, (seg) => slug(seg as Segment))
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("x1", (_d) => props.options.marginLeft)
            .attr("x2", (_d) => props.options.width - props.options.marginRight)
            .attr("y1", (d) => y2(d))
            .attr("y2", (d) => y2(d))
            .attr("stroke", (_d) => props.options.guideColour)
            .attr("stroke-width", 1),
        (update) => update,
        (exit) => exit.remove()
      );
  });
});
</script>

<template>
  <g>
    <g :x="group.winMin" :y="options.marginTop" :width="group.winWidth" :height="panelHeight" :fill="options.groupBackground" />
    <g :transform="axisTransform" ref="xAxis" />
    <g>
      <clipPath :id="boundingBoxId">
        <rect :x="group.winMin" :y="0" :width="group.winMax - group.winMin" :height="options.height"></rect>
      </clipPath>

      <g ref="guidesV1" :clip-path="'url(#' + boundingBoxId + ')'" />
      <g ref="guidesV2" :clip-path="'url(#' + boundingBoxId + ')'" />
      <g ref="guidesH1" :clip-path="'url(#' + boundingBoxId + ')'" />
      <g ref="guidesH2" :clip-path="'url(#' + boundingBoxId + ')'" />
    </g>
  </g>
</template>
