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

const boundingBoxId = computed<string>(() => {
  return `box-${props.group.winMin}-${props.group.winMax}-${props.options.height}`
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
      .data(props.group.segments, (seg) => slug(seg as Segment))
      .enter()
      .append("line")
      .attr("x1", (d) => x1(d))
      .attr("x2", (d) => x1(d))
      .attr("y1", (_d) => props.yScale(props.readMin))
      .attr("y2", (_d) => props.yScale(props.readMax))
      .attr("stroke", (_d) => props.options.guideColour)
      .attr("stroke-width", 1);
    d3.select(guidesV2.value)
      .selectAll("line")
      .data(props.group.segments, (seg) => slug(seg as Segment))
      .enter()
      .append("line")
      .attr("x1", (d) => x2(d))
      .attr("x2", (d) => x2(d))
      .attr("y1", (_d) => props.yScale(props.readMin))
      .attr("y2", (_d) => props.yScale(props.readMax))
      .attr("stroke", (_d) => props.options.guideColour)
      .attr("stroke-width", 1);

    d3.select(guidesH1.value)
      .selectAll("line")
      .data(props.group.segments, (seg) => slug(seg as Segment))
      .enter()
      .append("line")
      .attr("x1", (_d) => props.options.marginLeft)
      .attr("x2", (_d) => props.options.width - props.options.marginRight)
      .attr("y1", (d) => y1(d))
      .attr("y2", (d) => y1(d))
      .attr("stroke", (_d) => props.options.guideColour)
      .attr("stroke-width", 1);
    d3.select(guidesH2.value)
      .selectAll("line")
      .data(props.group.segments, (seg) => slug(seg as Segment))
      .enter()
      .append("line")
      .attr("x1", (_d) => props.options.marginLeft)
      .attr("x2", (_d) => props.options.width - props.options.marginRight)
      .attr("y1", (d) => y2(d))
      .attr("y2", (d) => y2(d))
      .attr("stroke", (_d) => props.options.guideColour)
      .attr("stroke-width", 1);
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

      <g ref="guidesV1" :clip-path="'url(#' + boundingBoxId + ')'"/>
      <g ref="guidesV2" :clip-path="'url(#' + boundingBoxId + ')'"/>
      <g ref="guidesH1" :clip-path="'url(#' + boundingBoxId + ')'"/>
      <g ref="guidesH2" :clip-path="'url(#' + boundingBoxId + ')'"/>
    </g>
  </g>
</template>
