<script setup lang="ts">
import { ReadInfo, Segment, SegmentGroupInfo } from "./segment";
import { Options } from "./options";
import { UnionFind } from "./unionfind";
import { computed, ComputedRef, onMounted, Ref, ref } from "vue";
import { computedAsync } from "@vueuse/core";
import ChromSegmentPlotBackground from "./ChromSegmentPlotBackground.vue";
import ChromSegmentPlotForeground from "./ChromSegmentPlotForeground.vue";
import * as d3 from "d3";
import { toPng } from "html-to-image";

const props = defineProps<{
  locus: string;
  chrom: string;
  segments: Segment[];
  options: Options;
  colours: Map<string, string>;
}>();

async function computeSha1(text: string): Promise<string> {
  const encoder = new TextEncoder();
  const data = encoder.encode(text);
  const hash = await window.crypto.subtle.digest("SHA-1", data);
  const hashArray = Array.from(new Uint8Array(hash)); // convert buffer to byte array
  const hashHex = hashArray.map((b) => b.toString(16).padStart(2, "0")).join(""); // convert bytes to hex string
  return hashHex;
}

const sig = computedAsync(() => {
  const readIds: string[] = d3.map(props.segments, (seg) => seg.readid);
  const digest = computeSha1(readIds.join());
  return digest;
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

    if (props.options.enableFlipping) {
      segs.forEach((seg) => {
        if (seg.pos == info.start && seg.offset > info.begin) {
          info.flip = true;
        }
      });
    }
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

  s.append("g")
    .attr("transform", `translate(12, ${props.options.height / 2})`)
    .append("text")
    .attr("text-anchor", "middle")
    .attr("transform", "rotate(-90)")
    .text("position in read");

  s.append("g")
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
    <v-btn icon="mdi-camera" class="float-right snapshot" @click="snap()"></v-btn>
    <svg :view-box="viewBox" :width="options.width" :height="options.height" ref="svg">
      <g ref="yAxis"></g>
      <ChromSegmentPlotBackground v-for="grp in segmentGroups" :key="sig + grp.grp"
      :y-scale="yScale"
      :group="grp"
      :read-info="readInfos"
      :read-min="readMin || 0"
      :read-max="readMax || 0"
      :options="options"
      :colours="colours"
      ></ChromSegmentPlotBackground>
      <ChromSegmentPlotForeground v-for="grp in segmentGroups" :key="sig + grp.grp"
      :y-scale="yScale"
      :group="grp"
      :read-info="readInfos"
      :read-min="readMin || 0"
      :read-max="readMax || 0"
      :options="options"
      :colours="colours"
      @selected-segment="(seg) => {clickedSegment = seg.readid; }"
      ></ChromSegmentPlotForeground>
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
