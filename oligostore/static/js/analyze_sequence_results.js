document.addEventListener("DOMContentLoaded", () => {
  document.querySelectorAll("[data-copy-target]").forEach((button) => {
    button.addEventListener("click", async () => {
      const targetId = button.dataset.copyTarget;
      const target = targetId ? document.getElementById(targetId) : null;
      if (!target) {
        return;
      }

      try {
        await navigator.clipboard.writeText(target.innerText);
      } catch (error) {
        window.alert("Failed to copy sequence to clipboard.");
      }
    });
  });
});
