document.addEventListener("DOMContentLoaded", () => {
  const button = document.getElementById("analyze-primerpair-btn");
  const config = document.getElementById("primerpair-combined-config");
  const resultDiv = document.getElementById("analysis-result");

  if (!button || !config || !resultDiv) {
    return;
  }

  button.addEventListener("click", async () => {
    const fwd = document.getElementById("id_forward_sequence")?.value || "";
    const rev = document.getElementById("id_reverse_sequence")?.value || "";
    const csrfToken = document.querySelector("[name=csrfmiddlewaretoken]")?.value || "";

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
        resultDiv.innerHTML = `<p class="text-error">${data.error}</p>`;
        return;
      }

      resultDiv.innerHTML = `
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
    } catch (error) {
      resultDiv.innerHTML = `<p class="text-error">${error.message}</p>`;
    }
  });
});
