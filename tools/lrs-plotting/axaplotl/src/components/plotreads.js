import * as d3 from "d3";

function sig(txt) {
  let hash = 50245149494948543n;
  for (let i = 0; i < txt.length; i += 1) {
    hash = (hash * 65978692084245217n + 237771813010709n * BigInt(txt.charCodeAt(i))) & 0xffffffffffffffffn;
  }
  let h0 = hash;
  let u = 0.0;
  for (let i = 0; i < 64; i += 1) {
    if ((hash & 1n) == 1n) {
      u += 1.0;
    }
    hash >>= 1n;
    u /= 2.0;
  }
  console.log(`${txt} -> ${h0} -> ${u}`);
  return u;
}

function plot_long_reads(rawData, elemId, options = {}) {
  console.log(rawData);
  
  const defaults = {
    width: 800,
    height: 600,
    marginTop: 50,
    marginRight: 50,
    marginBottom: 80,
    marginLeft: 80,
    gap: 30,
    groupBackground: "#fafafa",
    guideColour: "#f0f0f0",
    singleSegments: false,
    numTicks: 10,
    enableFlipping: true,
    randomizeColours: false,
  };
  options = { ...defaults, ...options };

  const root = d3.select(elemId);

  let dataByChrom = d3.group(rawData, (d) => d.chrom);
  dataByChrom.forEach((data, chrom) => {
    // First, group data by Read ID to figure out
    // whether to flip each one, and assign colours.
    //
    const groups = d3.group(data, (d) => d.readid);
    const groupInfo = new Map();
    let colourNumber = 0;
    groups.forEach((values, key) => {
      if (values.length == 1) {
        return;
      }
      const info = {
        start: d3.min(values, (d) => d.pos),
        length: d3.max(values, (d) => d.offset + d.qlen),
        begin: d3.min(values, (d) => d.offset),
        colourNumber: colourNumber,
        flip: false,
      };

      if (options.enableFlipping) {
        values.forEach((d) => {
          if (d.pos == info.start && d.offset > info.begin) {
            info.flip = true;
          }
        });
      }

      groupInfo.set(key, info);
      colourNumber += 1;
    });
    for (const info of groupInfo) {
      info[1].colourNumber /= colourNumber;
    }

    // Drop any reads that only have 1 segment.
    if (!options.singleSegments) {
      data = d3.filter(data, (d) => groupInfo.has(d.readid));
    }
    // Drop any segments that have zero rlen or dlen
    data = d3.filter(data, (d) => d.rlen > 0 && d.qlen > 0);

    if (data.length == 0) {
      return;
    }

    // The next phase is to group overlapping segments,
    // so we can identify chromosomal ranges that we want
    // to include on the x-axis. This is motivated by
    // the fact that sometimes we get large gaps, which
    // compress the axis, rendering the plots useless.
    // So we find overlaps, and drop the big gaps.
    //
    const uf = new UnionFind();
    const idx = new Map();
    for (let i = 0; i < data.length; i += 1) {
      const x = data[i];
      const xbegin = x.pos - 0.5 * x.rlen;
      const xend = x.pos + 1.5 * x.rlen;
      const xk = `x${i}`;
      uf.find(xk);
      idx.set(xk, x);

      for (let j = i + 1; j < data.length; j += 1) {
        const y = data[j];
        if (x.chrom != y.chrom) {
          // should never happen, since we already split by chrom
          continue;
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
    const segmentGroups = new Map();
    idx.forEach((x, xk) => {
      const gk = uf.find(xk);
      if (!segmentGroups.has(gk)) {
        segmentGroups.set(gk, [x]);
      } else {
        segmentGroups.get(gk).push(x);
      }
    });

    // Compile a list of segment positions for the groups.
    const segmentData = [];
    let totalWidth = 0;
    segmentGroups.forEach((items, grp) => {
      const grpMin = d3.min(items, (d) => d.pos);
      const grpMax = d3.max(items, (d) => d.pos + d.rlen);
      const grpWidth = grpMax - grpMin;
      totalWidth += grpWidth;
      segmentData.push({ grp, grpMin, grpMax, grpWidth });
    });
    segmentData.sort((a, b) => a.grpMin - b.grpMin);

    // Calculate how much space we have once we take out the gaps.
    const availableWidth = options.width - (options.marginLeft + options.marginRight + (segmentData.length - 1) * options.gap);

    // Work out the positions in the plot for each group.
    let leftOffset = options.marginLeft;
    segmentData.forEach((grp, i) => {
      grp.ratio = grp.grpWidth / totalWidth;
      grp.winMin = leftOffset;
      grp.winMax = leftOffset + grp.ratio * availableWidth;
      grp.winWidth = grp.ratio * availableWidth;
      grp.ticks = 1 + Math.floor(grp.ratio * options.numTicks);
      leftOffset += grp.ratio * availableWidth + options.gap;
    });

    const readMin = d3.min(data, (d) =>
      groupInfo.get(d.readid).flip ? groupInfo.get(d.readid).length - (d.offset + d.qlen) : d.offset
    );
    const readMax = d3.max(data, (d) =>
      groupInfo.get(d.readid).flip ? groupInfo.get(d.readid).length - d.offset : d.offset + d.qlen
    );
    // Declare the y (vertical position) scale.
    const yScale = d3
      .scaleLinear()
      .domain([readMin, readMax])
      .range([options.height - options.marginBottom, options.marginTop]);

    // Create the SVG container.
    const svg = d3.create("svg").attr("width", options.width).attr("height", options.height);

    const dotSize = 4;
    const defs = svg.append("svg:defs");
    const marker1 = defs
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

    const marker2 = defs
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

    // Add the title
    svg
      .append("g")
      .append("text")
      .text(chrom)
      .attr("transform", `translate(${options.marginLeft + 10}, ${options.marginTop - 10})`);

    function startMark(d) {
      const g = groupInfo.get(d.readid);
      if (d.offset == g.begin) {
        return "url(#dotClosed)";
      } else {
        return "";
      }
    }

    function endMark(d) {
      const g = groupInfo.get(d.readid);
      if (d.offset + d.qlen == g.length) {
        return "url(#dotOpen)";
      } else {
        return "";
      }
    }

    function y1(d) {
      let y = 0;
      if (d.strand == "+") {
        y = d.offset;
      } else {
        y = d.offset + d.qlen;
      }
      const info = groupInfo.get(d.readid);
      if (info.flip) {
        y = info.length - y;
      }
      return yScale(y);
    }

    function y2(d) {
      let y = 0;
      if (d.strand == "+") {
        y = d.offset + d.qlen;
      } else {
        y = d.offset;
      }
      const info = groupInfo.get(d.readid);
      if (info.flip) {
        y = info.length - y;
      }
      return yScale(y);
    }

    function findColour(d) {
      if (options.randomizeColours) {
        return d3.interpolateTurbo(sig(d.readid));
      }
      return d3.interpolateTurbo(groupInfo.get(d.readid).colourNumber);
    }

    // Add all the "background" elements
    //
    segmentData.forEach((grp, i) => {
      // Declare the x (horizontal position) scale.
      const xScale = d3.scaleLinear().domain([grp.grpMin, grp.grpMax]).range([grp.winMin, grp.winMax]);

      // Add a background
      svg
        .append("rect")
        .attr("x", grp.winMin)
        .attr("y", options.marginTop)
        .attr("width", grp.winWidth)
        .attr("height", options.height - (options.marginBottom + options.marginTop))
        .attr("fill", options.groupBackground);

      // Add the x-axis.
      svg
        .append("g")
        .attr("transform", `translate(0,${options.height - options.marginBottom})`)
        .call(d3.axisBottom(xScale).ticks(grp.ticks))
        .selectAll("text")
        .attr("y", 0)
        .attr("x", 9)
        .attr("dy", ".35em")
        .attr("transform", "rotate(45)")
        .style("text-anchor", "start");

      function x1(d) {
        return xScale(d.pos);
      }

      function x2(d) {
        return xScale(d.pos + d.rlen);
      }

      function dashes(d) {
        return "0";
      }

      const items = segmentGroups.get(grp.grp);

      // Add the vertical guilds
      svg
        .append("g")
        .selectAll("line")
        .data(items)
        .enter()
        .append("line")
        .attr("x1", (d) => x1(d))
        .attr("x2", (d) => x1(d))
        .attr("y1", (d) => yScale(readMin))
        .attr("y2", (d) => yScale(readMax))
        .attr("stroke", (d) => options.guideColour)
        .attr("stroke-width", 1);
      svg
        .append("g")
        .selectAll("line")
        .data(items)
        .enter()
        .append("line")
        .attr("x1", (d) => x2(d))
        .attr("x2", (d) => x2(d))
        .attr("y1", (d) => yScale(readMin))
        .attr("y2", (d) => yScale(readMax))
        .attr("stroke", (d) => options.guideColour)
        .attr("stroke-width", 1);

      // Add the horizontal guilds
      svg
        .append("g")
        .selectAll("line")
        .data(items)
        .enter()
        .append("line")
        .attr("x1", (d) => options.marginLeft)
        .attr("x2", (d) => options.width - options.marginRight)
        .attr("y1", (d) => y1(d))
        .attr("y2", (d) => y1(d))
        .attr("stroke", (d) => options.guideColour)
        .attr("stroke-width", 1);
      svg
        .append("g")
        .selectAll("line")
        .data(items)
        .enter()
        .append("line")
        .attr("x1", (d) => options.marginLeft)
        .attr("x2", (d) => options.width - options.marginRight)
        .attr("y1", (d) => y2(d))
        .attr("y2", (d) => y2(d))
        .attr("stroke", (d) => options.guideColour)
        .attr("stroke-width", 1);
    });

    // Add the actual segments.
    //
    segmentData.forEach((grp, i) => {
      // Declare the x (horizontal position) scale.
      const xScale = d3.scaleLinear().domain([grp.grpMin, grp.grpMax]).range([grp.winMin, grp.winMax]);

      function x1(d) {
        return xScale(+d.pos);
      }

      function x2(d) {
        return xScale(+d.pos + +d.rlen);
      }

      function dashes(d) {
        return "0";
      }

      const items = segmentGroups.get(grp.grp);

      // Add the lines
      svg
        .append("g")
        .selectAll("line")
        .data(items)
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

    // Add the y-axis.
    svg.append("g").attr("transform", `translate(${options.marginLeft},0)`).call(d3.axisLeft(yScale));

    // Append the SVG element.
    root.append(() => svg.node());
  });
}

class UnionFind {
  constructor() {
    this.parent = {};
    this.rank = {};
  }
  find(s) {
    if (!(s in this.parent)) {
      this.parent[s] = s;
      this.rank[s] = 0;
      return s;
    }
    const sp = this.parent[s];
    if (s !== sp) {
      this.parent[s] = this.find(sp);
    }
    return this.parent[s];
  }
  union(x, y) {
    const xr = this.find(x);
    const yr = this.find(y);
    if (xr === yr) {
      return xr;
    }
    let res = null;
    if (this.rank[xr] < this.rank[yr]) {
      this.parent[xr] = yr;
      res = yr;
    } else if (this.rank[xr] > this.rank[yr]) {
      this.parent[yr] = xr;
      res = xr;
    } else {
      this.parent[yr] = xr;
      this.rank[xr] += 1;
      res = xr;
    }
    return res;
  }
}

export default {
  plot_long_reads
}