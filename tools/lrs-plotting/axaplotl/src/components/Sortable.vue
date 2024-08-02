<script setup lang="ts" generic="T">
import { ref } from "vue";
import { ColumnSpecification, sortInPlace } from "./sortable";

defineProps<{
  columns: ColumnSpecification<T>[];
  keyGen: (item: T) => string;
  height: string;
}>();

const data = defineModel<T[]>("data");

const currentKey = ref<[string, boolean]>(["", true]);

function clickKey(spec: ColumnSpecification<T>) {
  if (currentKey.value) {
    if (currentKey.value[0] == spec.field.toString()) {
      if (currentKey.value[1] == true) {
        currentKey.value[1] = false;
      } else {
        return (currentKey.value = ["", true]);
      }
    } else {
      currentKey.value = [spec.field.toString(), true];
    }
  } else {
    currentKey.value = [spec.field.toString(), true];
  }
  sortInPlace(data.value || [], spec, currentKey.value[1]);
}
</script>

<template>
  <v-table fixed-header :height="height">
    <thead>
      <tr>
        <th v-for="spec in columns" :key="spec.field" :class="spec.kind == 'str' ? 'text-center' : 'text-right'" @click="clickKey(spec)">
          {{ spec.title }}
          <v-icon size="small" v-if="currentKey[0] == spec.field.toString()">{{
            currentKey[1] ? "mdi-chevron-up" : "mdi-chevron-down"
          }}</v-icon>
        </th>
      </tr>
    </thead>
    <tbody>
      <tr v-for="item in data" :key="keyGen(item)">
        <td
          v-for="spec in columns"
          :key="keyGen(item) + ':' + spec.field.toString()"
          :class="spec.kind == 'str' ? 'text-center' : 'text-right'"
        >
          {{ item[spec.field] }}
        </td>
      </tr>
    </tbody>
  </v-table>
</template>
