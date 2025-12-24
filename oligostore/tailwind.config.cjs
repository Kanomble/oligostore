import daisyui from "daisyui";

export default {
  content: [
    "./templates/**/*.html",
    "./core/templates/**/*.html",
    "./**/templates/**/*.html",
    "./**/*.py",
  ],

  theme: {
    extend: {},
  },

  plugins: [
    daisyui
  ],

  daisyui: {
    themes: [
      {
        light: {
          ...require("daisyui/src/theming/themes")["light"],
          "accent": "#3ABFF8",
          "accent-content": "#ffffff",
        }
      },
      "dark"
    ],
    logs: true,
  },
};
