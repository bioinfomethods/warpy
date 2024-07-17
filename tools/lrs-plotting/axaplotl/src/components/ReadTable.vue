<script setup lang="ts">
import { Ref } from "vue";
import { ReadItem } from "./segment";

defineProps<{
  items: ReadItem[];
}>();

const selected: Ref<string[]> = defineModel<string[]>({ required: true });

const headers = [
  { title: "Colour", value: "colour", align: "center", maxWidth: "2em" },
  { title: "Read ID", value: "readid", align: "start" },
];
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
    </v-data-table>
  </v-sheet>
</template>
