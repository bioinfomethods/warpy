<script setup lang="ts">
import { Segment } from "./segment";
import { Options } from "./options";
import { UnionFind } from "./unionfind";
import { computed, ComputedRef, onMounted, ref } from "vue";
import * as d3 from "d3";
import { toPng } from "html-to-image";

const props = defineProps<{
  chrom: string;
  segments: Segment[];
  options: Options;
  colours: Map<string, number>;
}>();

type ReadInfo = {
  start: number;
  length: number;
  begin: number;
  flip: boolean;
};

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

type SegmentGroupInfo = {
  grp: string;
  grpMin: number;
  grpMax: number;
  grpWidth: number;
  ratio: number;
  winMin: number;
  winMax: number;
  winWidth: number;
  ticks: number;
  segments: Segment[];
};

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

const xScales = computed(() => {
  const res = segmentGroups.value.map((grp) => {
    return d3.scaleLinear().domain([grp.grpMin, grp.grpMax]).range([grp.winMin, grp.winMax]);
  });
  return res;
});

function findColour(seg: Segment) {
  const col = props.colours.get(seg.readid) || 0;
  return d3.interpolateTurbo(col);
}

function startMark(seg: Segment) {
  const g = readInfos.value.get(seg.readid);
  if (g && seg.offset == g.begin) {
    return "url(#dotClosed)";
  } else {
    return "";
  }
}

function endMark(seg: Segment): string {
  const g = readInfos.value.get(seg.readid);
  if (g && seg.offset + seg.qlen == g.length) {
    return "url(#dotOpen)";
  } else {
    return "";
  }
}

const svg = ref(null);
const yAxis = ref(null);
const groups = ref([]);

const viewBox = computed(() => {
  const w = props.options.width;
  const h = props.options.height;
  return `0 0 ${w} ${h}`;
});

onMounted(() => {
  const opts = props.options;

  const s = d3.select(svg.value);

  const dotSize = 4;
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

  if (false) {
    s.append("g")
      .append("text")
      .text(props.chrom)
      .attr("transform", `translate(${props.options.marginLeft + 10}, ${props.options.marginTop - 10})`);
  }

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

  const gy = d3.select(yAxis.value);
  const leftAxis = d3.axisLeft(yScale.value);
  gy.attr("transform", `translate(${props.options.marginLeft}, 0)`).call(leftAxis as any);

  function y1(seg: Segment): number {
    let y = 0;
    if (seg.strand == "+") {
      y = seg.offset;
    } else {
      y = seg.offset + seg.qlen;
    }
    const info = readInfos.value.get(seg.readid);
    if (info && info.flip) {
      y = info.length - y;
    }
    return yScale.value(y);
  }

  function y2(seg: Segment): number {
    let y = 0;
    if (seg.strand == "+") {
      y = seg.offset + seg.qlen;
    } else {
      y = seg.offset;
    }
    const info = readInfos.value.get(seg.readid);
    if (info && info.flip) {
      y = info.length - y;
    }
    return yScale.value(y);
  }

  groups.value.forEach((e, i) => {
    const grp = segmentGroups.value[i];
    /*
     */
    const xScale = xScales.value[i];

    function x1(seg: Segment): number {
      return xScale(seg.pos);
    }

    function x2(seg: Segment) {
      return xScale(seg.pos + seg.rlen);
    }

    function dashes(_seg: Segment): string {
      return "0";
    }

    const gg = d3.select(e);
    gg.append("rect")
      .attr("x", grp.winMin)
      .attr("y", opts.marginTop)
      .attr("width", grp.winWidth)
      .attr("height", opts.height - (opts.marginBottom + opts.marginTop))
      .attr("fill", opts.groupBackground);

    gg.append("g")
      .attr("transform", `translate(0,${opts.height - opts.marginBottom})`)
      .call(d3.axisBottom(xScale).ticks(grp.ticks))
      .selectAll("text")
      .attr("y", 0)
      .attr("x", 9)
      .attr("dy", ".35em")
      .attr("transform", "rotate(45)")
      .style("text-anchor", "start");

    gg.append("g")
      .selectAll("line")
      .data(grp.segments)
      .enter()
      .append("line")
      .attr("x1", (d) => x1(d))
      .attr("x2", (d) => x1(d))
      .attr("y1", (_d) => yScale.value(readMin.value || 0))
      .attr("y2", (_d) => yScale.value(readMax.value || 0))
      .attr("stroke", (_d) => opts.guideColour)
      .attr("stroke-width", 1);
    gg.append("g")
      .selectAll("line")
      .data(grp.segments)
      .enter()
      .append("line")
      .attr("x1", (d) => x2(d))
      .attr("x2", (d) => x2(d))
      .attr("y1", (_d) => yScale.value(readMin.value || 0))
      .attr("y2", (_d) => yScale.value(readMax.value || 0))
      .attr("stroke", (_d) => opts.guideColour)
      .attr("stroke-width", 1);

    gg.append("g")
      .selectAll("line")
      .data(grp.segments)
      .enter()
      .append("line")
      .attr("x1", (_d) => opts.marginLeft)
      .attr("x2", (_d) => opts.width - opts.marginRight)
      .attr("y1", (d) => y1(d))
      .attr("y2", (d) => y1(d))
      .attr("stroke", (_d) => opts.guideColour)
      .attr("stroke-width", 1);
    gg.append("g")
      .selectAll("line")
      .data(grp.segments)
      .enter()
      .append("line")
      .attr("x1", (_d) => opts.marginLeft)
      .attr("x2", (_d) => opts.width - opts.marginRight)
      .attr("y1", (d) => y2(d))
      .attr("y2", (d) => y2(d))
      .attr("stroke", (_d) => opts.guideColour)
      .attr("stroke-width", 1);

    gg.append("g")
      .selectAll("line")
      .data(grp.segments)
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
      .style("stroke-dasharray", (d) => dashes(d));
  });
});

async function snap() {
  const svgElem = svg.value;
  if (svgElem) {
    const dataUrl = await toPng(svgElem);
    const img = new Image();
    img.src = dataUrl;
    //download(dataUrl, 'plot.png');
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
      <g v-for="grp in segmentGroups" :key="grp.grp" ref="groups"></g>
    </svg>
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
