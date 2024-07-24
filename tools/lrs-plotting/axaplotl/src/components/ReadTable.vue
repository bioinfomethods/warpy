<script setup lang="ts">
import { Ref } from "vue";
import { Locus, ReadItem } from "./segment";

defineProps<{}>();

const selected: Ref<string[]> = defineModel<string[]>("selected", { required: true });
const items: Ref<ReadItem[]> = defineModel<ReadItem[]>("items", { required: true });

const headers = [
  { title: "Colour", value: "colour", align: "center", maxWidth: "4em" },
  { title: "Flipped", value: "flipped", align: "center", maxWidth: "4em" },
  { title: "Strand", value: "strand", align: "center", maxWidth: "4em" },
  { title: "Length", value: "length", align: "end"}, 
  { title: "Read ID", value: "readid", align: "center" },
  { title: "Mapped Locations", value: "mapped", align: "start" },
];

function formatLocus(locus: Locus): string {
  return `${locus.chrom}:${locus.start}-${locus.end}`;
}
</script>

<template>
  <v-sheet rounded border>
    <v-data-table
      :items="items"
      :headers="headers as any"
      density="compact"
      items-per-page="-1"
      hide-default-footer
      show-select
      v-model="selected"
      item-value="readid"
      item-key="readid"
    >
      <template v-slot:item.colour="{ item }">
        <v-icon icon="mdi-circle" :color="item.colour"></v-icon>
      </template>
      <template v-slot:item.flipped="{ item }">
        <v-checkbox-btn density="compact" v-model="item.flipped"></v-checkbox-btn>
      </template>
      <template v-slot:item.mapped="{ item }">
        <div class="text-left text-caption">
          <span v-for="(locus, i) in item.mapped">
            <span v-if="i > 0">, </span>
            <span class="nowrap">
              {{ formatLocus(locus) }}
            </span>
          </span>
        </div>
      </template>
    </v-data-table>
  </v-sheet>
</template>

<style scoped>
.leftify {
  width: 100%;
}
.nowrap {
  white-space: nowrap;
}
</style>
