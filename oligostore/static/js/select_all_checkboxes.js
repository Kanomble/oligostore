document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("[data-select-all-target]").forEach((toggle) => {
    const selector = toggle.dataset.selectAllTarget;
    if (!selector) {
      return;
    }

    const syncItems = () => {
      document.querySelectorAll(selector).forEach((checkbox) => {
        checkbox.checked = toggle.checked;
      });
    };

    toggle.addEventListener("change", syncItems);
  });
});
