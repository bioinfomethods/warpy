<script setup lang="ts">
import { Locus, parseLocus, ReadInfo, ReadItem, Segment, SegmentGroupInfo, validLocus } from "./segment";
import { Options } from "./options";
import { UnionFind } from "./unionfind";
import { computeSha1, humanize } from "./utils";
import { computed, ComputedRef, onMounted, Ref, ref, watchEffect } from "vue";
import { computedAsync } from "@vueuse/core";
import ChromSegmentPlotBackground from "./ChromSegmentPlotBackground.vue";
import ChromSegmentPlotForeground from "./ChromSegmentPlotForeground.vue";
import * as d3 from "d3";
import { toPng } from "html-to-image";

const props = defineProps<{
  locus: string;
  chrom: string;
  segments: Segment[];
  reads: ReadItem[];
  options: Options;
  colours: Map<string, string>;
}>();

const scaleAutoSwitch = ref<boolean>(true);

const sig = computedAsync(() => {
  const readIds: string[] = d3.map(props.segments, (seg) => seg.readid);
  const digest = computeSha1(readIds.join());
  return digest;
});

const refMin = computed<number>(() => {
  return d3.min(props.segments, (seg) => seg.pos) || 0;
});

const refMax = computed<number>(() => {
  return d3.max(props.segments, (seg) => seg.pos + seg.rlen) || 0;
});

const refRange = computed<number>(() => {
  if (refMin.value && refMax.value) {
    return refMax.value - refMin.value;
  }
  return 0;
});

const dataLocusText = computed<string>(() => {
  if (refMin.value != 0 && refMax.value != 0) {
    return `${props.chrom}:${refMin.value}-${refMax.value}`;
  } else {
    return "";
  }
});

const dataLocus = computed<Locus | null>(() => {
  if (dataLocusText.value) {
    return parseLocus(dataLocusText.value);
  }
  return null;
});

const locusBoxText = ref<string>("");

const locusText = ref<string>("");

watchEffect(() => {
  if (scaleAutoSwitch.value == false && (locusBoxText.value == undefined || locusBoxText.value == "") && dataLocusText.value) {
    locusBoxText.value = dataLocusText.value;
  }

  if (locusBoxText.value && validLocus(locusBoxText.value) == true) {
    locusText.value = locusBoxText.value;
  }
});

const locusValue = computed<Locus | null>(() => {
  if (locusText.value) {
    return parseLocus(locusText.value);
  }
  return null;
});

const zoomOutScale = 1.1;

function zoomIn() {
  if (locusValue.value != null) {
    const chrom = locusValue.value.chrom;
    let start = locusValue.value.start;
    let end = locusValue.value.end;
    const middle = (start + end) / 2;
    let range = end - start;
    range /= zoomOutScale;
    if (range < 100 || (refRange.value && range > 1.05 * refRange.value)) {
      return;
    }
    start = Math.floor(middle - 0.5 * range);
    end = Math.ceil(middle + 0.5 * range);
    const locus = `${chrom}:${start}-${end}`;
    locusBoxText.value = locus;
  }
}

function zoomOut() {
  if (locusValue.value != null) {
    const chrom = locusValue.value.chrom;
    let start = locusValue.value.start;
    let end = locusValue.value.end;
    const middle = (start + end) / 2;
    let range = end - start;
    range *= zoomOutScale;
    if (range < 100 || (refRange.value && range > 1.05 * refRange.value)) {
      return;
    }
    start = Math.floor(middle - 0.5 * range);
    end = Math.ceil(middle + 0.5 * range);
    const locus = `${chrom}:${start}-${end}`;
    locusBoxText.value = locus;
  }
}

const panScale = 0.1;

function panLeft() {
  if (locusValue.value != null) {
    const chrom = locusValue.value.chrom;
    let start = locusValue.value.start;
    let end = locusValue.value.end;
    let range = end - start;
    const shift = Math.round(range * panScale);
    if (shift == 0) {
      return;
    }
    start -= shift;
    end -= shift;
    if (start < 0) {
      return;
    }
    const locus = `${chrom}:${start}-${end}`;
    locusBoxText.value = locus;
  }
}

function panRight() {
  if (locusValue.value != null) {
    const chrom = locusValue.value.chrom;
    let start = locusValue.value.start;
    let end = locusValue.value.end;
    let range = end - start;
    const shift = Math.round(range * panScale);
    if (shift == 0) {
      return;
    }
    start += shift;
    end += shift;
    const locus = `${chrom}:${start}-${end}`;
    locusBoxText.value = locus;
  }
}

const everythingSegmentGroup = computed<SegmentGroupInfo>(() => {
  const locus = scaleAutoSwitch.value ? dataLocus.value : locusValue.value;
  const grpMin: number = locus?.start || 0;
  const grpMax: number = locus?.end || 0;
  const grpWidth = grpMax - grpMin;
  const res: SegmentGroupInfo = {
    grp: "all",
    grpMin,
    grpMax,
    grpWidth,
    ratio: 1.0,
    winMin: props.options.marginLeft,
    winMax: props.options.width - props.options.marginRight,
    winWidth: props.options.marginRight - props.options.marginLeft,
    ticks: props.options.numTicks,
    segments: props.segments,
  };
  return res;
});

const readInfos = computed(() => {
  const res: Map<string, ReadInfo> = new Map();
  const groups = d3.group(props.segments, (seg) => seg.readid);
  groups.forEach((segs, readid) => {
    const info: ReadInfo = {
      start: d3.min(segs, (d) => d.pos) || 0,
      length: d3.max(segs, (d) => d.offset + d.qlen) || 0,
      begin: d3.min(segs, (d) => d.offset) || 0,
      flip: false,
    };
    res.set(readid, info);
  });
  return res;
});

function getReadInfo(readid: string): ReadInfo {
  const r = readInfos.value;
  const ri = r.get(readid);
  if (ri) {
    return ri;
  } else {
    throw `readid: ${readid} not found`;
  }
}

const readIndex = computed<Map<string,ReadItem>>(() => {
  const res: Map<string,ReadItem> = new Map();
  props.reads.forEach((item) => {
    res.set(item.readid, item);
  });
  return res;
});


const segmentGroups = computed(() => {
  const uf = new UnionFind();
  const idx: Map<string, Segment> = new Map();

  for (let i = 0; i < props.segments.length; i += 1) {
    const x = props.segments[i];
    const xbegin = x.pos - 0.5 * x.rlen;
    const xend = x.pos + 1.5 * x.rlen;
    const xk = `x${i}`;
    uf.find(xk);
    idx.set(xk, x);

    for (let j = i + 1; j < props.segments.length; j += 1) {
      const y = props.segments[j];
      if (x.chrom != y.chrom) {
        throw `paritioning by chromosome failed! ${x.chrom} ${y.chrom}`;
      }
      const ybegin = y.pos - 0.5 * y.rlen;
      const yend = y.pos + 1.5 * y.rlen;
      const yk = `x${j}`;
      uf.find(yk);

      if (xbegin <= yend && ybegin <= xend) {
        uf.union(xk, yk);
      }
    }
  }

  const segmentGroups: Map<string, Segment[]> = new Map();
  idx.forEach((x, xk) => {
    const gk = uf.find(xk);
    const grp = segmentGroups.get(gk);
    if (grp) {
      grp.push(x);
    } else {
      segmentGroups.set(gk, [x]);
    }
  });

  const res: SegmentGroupInfo[] = [];
  let totalWidth = 0;
  segmentGroups.forEach((group, gk) => {
    const grpMin = d3.min(group, (d) => d.pos) || 0;
    const grpMax = d3.max(group, (d) => d.pos + d.rlen) || 0;
    const grpWidth = grpMax - grpMin;
    totalWidth += grpWidth;
    res.push({
      grp: gk,
      grpMin,
      grpMax,
      grpWidth,
      ratio: 0,
      winMin: 0,
      winMax: 0,
      winWidth: 0,
      ticks: 0,
      segments: group,
    });
  });
  res.sort((a, b) => a.grpMin - b.grpMin);

  const opts = props.options;
  const availableWidth = opts.width - (opts.marginLeft + opts.marginRight + (res.length - 1) * opts.gap);
  let leftOffset = opts.marginLeft;
  res.forEach((grp, _i) => {
    grp.ratio = grp.grpWidth / totalWidth;
    grp.winMin = leftOffset;
    grp.winMax = leftOffset + grp.ratio * availableWidth;
    grp.winWidth = grp.ratio * availableWidth;
    grp.ticks = 1 + Math.floor(grp.ratio * opts.numTicks);
    leftOffset += grp.ratio * availableWidth + opts.gap;
  });

  return res;
});

const readMin = computed(() => {
  return d3.min(props.segments, (seg) =>
    getReadInfo(seg.readid).flip ? getReadInfo(seg.readid).length - (seg.offset + seg.qlen) : seg.offset
  );
});

const readMax = computed(() => {
  return d3.max(props.segments, (seg) =>
    getReadInfo(seg.readid).flip ? getReadInfo(seg.readid).length - seg.offset : seg.offset + seg.qlen
  );
});

const yScale: ComputedRef<d3.ScaleLinear<number, number, never>> = computed(() => {
  const rLo = readMin.value || 0;
  const rHi = readMax.value || 0;
  const opts = props.options;
  return d3
    .scaleLinear()
    .domain([rLo, rHi])
    .range([opts.height - opts.marginBottom, opts.marginTop]);
});

const svg = ref(null);
const yAxis = ref(null);
const yAxisTitle = ref(null);
const xAxisTitle = ref(null);

const viewBox = computed(() => {
  const w = props.options.width;
  const h = props.options.height;
  return `0 0 ${w} ${h}`;
});

const clickedSegment: Ref<string> = ref("");

async function copyToClipboard() {
  await navigator.clipboard.writeText(clickedSegment.value);
}

onMounted(() => {
  const s = d3.select(svg.value);

  const dotSize = 5;
  const defs = s.append("svg:defs");
  defs
    .append("svg:marker")
    .attr("id", "dotClosed")
    .attr("viewBox", [0, 0, 20, 20])
    .attr("refX", dotSize)
    .attr("refY", dotSize)
    .attr("markerWidth", dotSize)
    .attr("markerHeight", dotSize)
    .append("circle")
    .attr("cx", dotSize)
    .attr("cy", dotSize)
    .attr("r", dotSize)
    .style("fill", "grey");

  defs
    .append("svg:marker")
    .attr("id", "dotOpen")
    .attr("viewBox", [0, 0, 20, 20])
    .attr("refX", dotSize)
    .attr("refY", dotSize)
    .attr("markerWidth", dotSize)
    .attr("markerHeight", dotSize)
    .append("circle")
    .attr("cx", dotSize)
    .attr("cy", dotSize)
    .attr("r", dotSize)
    .attr("stroke-width", 1)
    .attr("stroke", "grey")
    .style("fill", "none");

  d3.select(yAxis.value).attr("transform", `translate(${props.options.marginLeft},0)`).call(d3.axisLeft(yScale.value) as any);

  d3.select(yAxisTitle.value)
    .attr("transform", `translate(12, ${props.options.height / 2})`)
    .append("text")
    .attr("text-anchor", "middle")
    .attr("transform", "rotate(-90)")
    .text("position in read");

  d3.select(xAxisTitle.value)
    .attr("transform", `translate(${props.options.width / 2}, ${props.options.height - 12})`)
    .append("text")
    .attr("text-anchor", "middle")
    .text(`position in reference (${props.chrom})`);
});

async function snap() {
  const svgElem = svg.value;
  if (svgElem) {
    const dataUrl = await toPng(svgElem);
    const img = new Image();
    img.src = dataUrl;
    const win = window.open();
    if (win) {
      win.document.body.style.width = "100%";
      win.document.body.style.height = "100%";
      win.document.body.appendChild(img);
    }
  }
}
</script>

<template>
  <div>
    <v-card>
      <v-card-actions>
        <v-switch label="Auto scale" v-model="scaleAutoSwitch" color="primary"></v-switch>
        <v-text-field
          label="Locus"
          :disabled="scaleAutoSwitch"
          density="compact"
          max-width="32em"
          v-model="locusBoxText"
          :rules="[validLocus]"
          prepend-inner-icon="mdi-creation"
          @click:prepend-inner="locusBoxText = dataLocusText"
          append-inner-icon="mdi-refresh"
        >
        </v-text-field>
        <v-btn icon="mdi-plus-circle-outline" density="compact" :disabled="scaleAutoSwitch" @click="zoomIn()"></v-btn>
        <v-btn icon="mdi-minus-circle-outline" density="compact" :disabled="scaleAutoSwitch" @click="zoomOut()"></v-btn>
        <v-btn icon="mdi-arrow-left-bold-circle-outline" density="compact" :disabled="scaleAutoSwitch" @click="panLeft()"></v-btn>
        <v-btn icon="mdi-arrow-right-bold-circle-outline" density="compact" :disabled="scaleAutoSwitch" @click="panRight()"></v-btn>
      </v-card-actions>
      <span v-if="locusValue != null">{{ humanize(locusValue.end - locusValue.start) }}</span>
    </v-card>
    <v-btn icon="mdi-camera" class="float-right snapshot" @click="snap()"></v-btn>
    <svg :view-box="viewBox" :width="options.width" :height="options.height" ref="svg">
      <g ref="yAxis"></g>
      <g ref="yAxisTitle"></g>
      <g ref="xAxisTitle"></g>
      <g v-if="scaleAutoSwitch">
        <ChromSegmentPlotBackground
          v-for="grp in segmentGroups"
          :key="sig + grp.grp"
          :y-scale="yScale"
          :group="grp"
          :reads="readIndex"
          :read-min="readMin || 0"
          :read-max="readMax || 0"
          :options="options"
          :colours="colours"
        ></ChromSegmentPlotBackground>
        <ChromSegmentPlotForeground
          v-for="grp in segmentGroups"
          :key="sig + grp.grp"
          :y-scale="yScale"
          :group="grp"
          :reads="readIndex"
          :read-min="readMin || 0"
          :read-max="readMax || 0"
          :options="options"
          :colours="colours"
          @selected-segment="
            (seg) => {
              clickedSegment = seg.readid;
            }
          "
        ></ChromSegmentPlotForeground>
      </g>
      <g v-else>
        <ChromSegmentPlotBackground
        :key="locusText"
          :y-scale="yScale"
          :group="everythingSegmentGroup"
          :reads="readIndex"
          :read-min="readMin || 0"
          :read-max="readMax || 0"
          :options="options"
          :colours="colours"
        ></ChromSegmentPlotBackground>
        <ChromSegmentPlotForeground
        :key="locusText"
          :y-scale="yScale"
          :group="everythingSegmentGroup"
          :reads="readIndex"
          :read-min="readMin || 0"
          :read-max="readMax || 0"
          :options="options"
          :colours="colours"
          @selected-segment="
            (seg) => {
              clickedSegment = seg.readid;
            }
          "
        ></ChromSegmentPlotForeground>
      </g>
    </svg>
    <h3>
      {{ clickedSegment }}
      <v-btn size="x-small" variant="plain" icon="mdi-clipboard-text" v-if="clickedSegment" @click="copyToClipboard"></v-btn>
    </h3>
  </div>
</template>

<style scoped>
svg {
  background-color: white;
}
.snapshot {
  margin: 5pt;
}
</style>
