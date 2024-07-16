export class UnionFind {
    parent: { [key: string]: string };
    rank: { [key: string]: number };
  
    constructor() {
      this.parent = {};
      this.rank = {};
    }
    find(s: string): string {
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
    union(x: string, y: string): string {
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