<script setup lang="ts">
import { flipSegmentIfNecessary, ReadItem, Segment, SegmentGroupInfo } from "./segment";
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

const verticalGuidePositions = computed<number[]>(() => {
  const res: number[] = [];
  if (segments.value) {
    segments.value.forEach((seg) => {
      res.push(x1(seg));
      res.push(x2(seg));
    });
  }
  return [...new Set(res)];
});

const horizontalGuidePositions = computed<number[]>(() => {
  const res: number[] = [];
  if (segments.value) {
    segments.value.forEach((seg) => {
      res.push(y1(seg));
      res.push(y2(seg));
    });
  }
  return [...new Set(res)];
});
const xAxis = ref(null);
const guidesV = ref(null);
const guidesH = ref(null);

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

    d3.select(guidesV.value)
      .selectAll("line")
      .data(verticalGuidePositions.value, (n) => (n as number).toString())
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("x1", (d) => d)
            .attr("x2", (d) => d)
            .attr("y1", (_d) => props.yScale(props.readMin))
            .attr("y2", (_d) => props.yScale(props.readMax))
            .attr("stroke", (_d) => props.options.guideColour)
            .attr("stroke-width", 1),
        (update) => update,
        (exit) => exit.remove()
      );

    d3.select(guidesH.value)
      .selectAll("line")
      .data(horizontalGuidePositions.value, (n) => (n as number).toString())
      .join(
        (enter) =>
          enter
            .append("line")
            .attr("x1", (_d) => props.options.marginLeft)
            .attr("x2", (_d) => props.options.width - props.options.marginRight)
            .attr("y1", (d) => d)
            .attr("y2", (d) => d)
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

      <g ref="guidesV" :clip-path="'url(#' + boundingBoxId + ')'" />
      <g ref="guidesH" :clip-path="'url(#' + boundingBoxId + ')'" />
    </g>
  </g>
</template>
