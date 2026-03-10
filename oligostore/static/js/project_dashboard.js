document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("[data-toggle-target]").forEach((button) => {
    button.addEventListener("click", () => {
      const target = document.getElementById(button.dataset.toggleTarget);
      if (!target) {
        return;
      }
      target.classList.toggle("hidden");
      button.setAttribute("aria-expanded", String(!target.classList.contains("hidden")));
    });
  });

  document.querySelectorAll("[data-row-link]").forEach((row) => {
    row.addEventListener("click", (event) => {
      if (event.target.closest("a, button")) {
        return;
      }
      const href = row.dataset.rowLink;
      if (!href) {
        return;
      }
      const confirmMessage = row.dataset.confirm;
      if (confirmMessage && !window.confirm(confirmMessage)) {
        return;
      }
      window.location.href = href;
    });
  });
});
