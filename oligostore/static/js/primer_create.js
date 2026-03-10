document.addEventListener("DOMContentLoaded", () => {
  const button = document.getElementById("analyze-primer-btn");
  const config = document.getElementById("primer-create-config");
  const resultDiv = document.getElementById("analysis-result");

  if (!button || !config || !resultDiv) {
    return;
  }

  button.addEventListener("click", async () => {
    const seq = document.getElementById("id_sequence")?.value || "";
    const csrfToken = document.querySelector("[name=csrfmiddlewaretoken]")?.value || "";

    try {
      const response = await fetch(config.dataset.analyzePrimerUrl, {
        method: "POST",
        headers: {
          "X-CSRFToken": csrfToken,
          "Content-Type": "application/x-www-form-urlencoded",
        },
        body: new URLSearchParams({ sequence: seq }),
      });
      const data = await response.json();

      if (data.error) {
        resultDiv.innerHTML = `<p class="text-error">${data.error}</p>`;
        return;
      }

      resultDiv.innerHTML = `
        <div class="card bg-base-200 shadow p-4">
          <h3 class="font-bold text-lg mb-2">Analysis Result</h3>
          <p class="strong">Tm: ${data.tm}</p>
          <p class="strong">GC: ${data.gc_content}</p>
          <p class="strong">Hairpin Found: ${data.hairpin}</p>
          <p class="strong">Hairpin: ${data.hairpin_dg}</p>
          <p class="strong">Self Dimer Found: ${data.self_dimer}</p>
          <p class="strong">Self-Dimer: ${data.self_dimer_dg}</p>
        </div>
      `;
    } catch (error) {
      resultDiv.innerHTML = `<p class="text-error">${error.message}</p>`;
    }
  });
});
