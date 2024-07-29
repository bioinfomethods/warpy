<script setup lang="ts">
import StaticPlot from "./components/StaticPlot.vue";
import DynamicPlot from "./components/DynamicPlot.vue";
import { Options } from "./components/options";

const options: Partial<Options> = {};

const staticOrDynamic = defineModel<string>("staticOrDynamic");
</script>

<template>
  <v-sheet>
    <v-tabs v-model="staticOrDynamic">
      <v-tab value="static">Precomputed Alignment Summaries</v-tab>
      <v-tab value="dynamic">Alignment Summaries from BAM</v-tab>
    </v-tabs>
    <v-tabs-window v-model="staticOrDynamic">
      <v-tabs-window-item value="static">
        <StaticPlot :options="options"></StaticPlot>
      </v-tabs-window-item>
      <v-tabs-window-item value="dynamic">
        <DynamicPlot :options="options"></DynamicPlot>
      </v-tabs-window-item>
    </v-tabs-window>
  </v-sheet>
</template>
