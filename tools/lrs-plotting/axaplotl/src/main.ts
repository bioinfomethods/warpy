import { createApp } from "vue";
import "vuetify/styles";
import { createVuetify } from "vuetify";
import * as components from "vuetify/components";
import * as directives from "vuetify/directives";

import "./style.css";
import '@mdi/font/css/materialdesignicons.css';

import App from "./App.vue";

const vuetify = createVuetify({
  ssr: true,
  components,
  directives,
});

createApp(App).use(vuetify).mount("#app");
