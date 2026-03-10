document.addEventListener("DOMContentLoaded", () => {
  const config = document.getElementById("primerpair-list-config");
  if (!config) {
    return;
  }

  const analysisButtons = document.querySelectorAll(".analyze-primerpair-btn");
  const csrfToken = document.querySelector("[name=csrfmiddlewaretoken]")?.value || "";

  analysisButtons.forEach((button) => {
    button.addEventListener("click", async () => {
      const fwd = button.dataset.forward || "";
      const rev = button.dataset.reverse || "";
      const targetId = button.dataset.target;
      const targetRow = targetId ? document.getElementById(targetId) : null;
      const content = targetRow?.querySelector(".analysis-content");

      if (!content || !csrfToken) {
        return;
      }

      targetRow.classList.remove("hidden");
      content.innerHTML = "<p class=\"text-sm opacity-70\">Analyzing primer compatibility...</p>";

      try {
        const response = await fetch(config.dataset.analyzePrimerpairUrl, {
          method: "POST",
          headers: {
            "X-CSRFToken": csrfToken,
            "Content-Type": "application/x-www-form-urlencoded",
          },
          body: new URLSearchParams({
            forward_sequence: fwd,
            reverse_sequence: rev,
          }),
        });
        const data = await response.json();

        if (data.error) {
          content.innerHTML = `<p class="text-error">${data.error}</p>`;
          return;
        }

        content.innerHTML = `
          <div class="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div class="card bg-base-200 p-4">
              <h3 class="font-bold text-lg mb-2">Forward Primer</h3>
              <p><strong>Tm:</strong> ${data.forward.tm}</p>
              <p><strong>GC:</strong> ${data.forward.gc}%</p>
              <p><strong>Hairpin dG:</strong> ${data.forward.hairpin_dg}</p>
              <p><strong>Self-Dimer dG:</strong> ${data.forward.self_dimer_dg}</p>
            </div>
            <div class="card bg-base-200 p-4">
              <h3 class="font-bold text-lg mb-2">Reverse Primer</h3>
              <p><strong>Tm:</strong> ${data.reverse.tm}</p>
              <p><strong>GC:</strong> ${data.reverse.gc}%</p>
              <p><strong>Hairpin dG:</strong> ${data.reverse.hairpin_dg}</p>
              <p><strong>Self-Dimer dG:</strong> ${data.reverse.self_dimer_dg}</p>
            </div>
          </div>
          <div class="card bg-base-300 p-4 mt-4">
            <h3 class="font-bold text-lg mb-2">Hetero-Dimer (Forward / Reverse)</h3>
            <p><strong>dG:</strong> ${data.hetero_dimer_dg}</p>
          </div>
        `;
      } catch {
        content.innerHTML = "<p class=\"text-error\">Analysis failed. Please try again.</p>";
      }
    });
  });
});
