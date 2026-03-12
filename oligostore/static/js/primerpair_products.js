document.addEventListener("DOMContentLoaded", () => {
  const form = document.getElementById("primerpair-products-form");
  const config = document.getElementById("primerpair-products-config");
  const statusBox = document.getElementById("primerpair-products-status");
  const errorBox = document.getElementById("primerpair-products-error");
  const resultsBox = document.getElementById("primerpair-products-results");
  const summaryBox = document.getElementById("primerpair-products-summary");
  const emptyBox = document.getElementById("primerpair-products-empty");
  const listBox = document.getElementById("primerpair-products-list");

  if (!form || !config || !statusBox || !errorBox || !resultsBox || !summaryBox || !emptyBox || !listBox) {
    return;
  }

  const errorMessage = errorBox.querySelector(".alert");
  const asyncUrl = config.dataset.asyncUrl;
  const statusUrlTemplate = config.dataset.statusUrlTemplate;
  const downloadUrl = config.dataset.downloadUrl;

  function getCookie(name) {
    const value = `; ${document.cookie}`;
    const parts = value.split(`; ${name}=`);
    if (parts.length === 2) {
      return parts.pop().split(";").shift();
    }
    return "";
  }

  function escapeHtml(value) {
    return String(value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  function resetUI() {
    statusBox.classList.add("hidden");
    errorBox.classList.add("hidden");
    resultsBox.classList.add("hidden");
    emptyBox.classList.add("hidden");
    summaryBox.textContent = "";
    listBox.innerHTML = "";
  }

  function showError(message) {
    errorMessage.textContent = message;
    errorBox.classList.remove("hidden");
  }

  function renderProducts(result) {
    const primerPair = result.primer_pair || {};
    const sequenceFile = result.sequence_file || {};
    const products = result.products || [];

    resultsBox.classList.remove("hidden");
    summaryBox.textContent = `Showing results for ${primerPair.name || "selected pair"} against ${sequenceFile.name || "selected sequence file"}.`;

    if (!products.length) {
      emptyBox.classList.remove("hidden");
      return;
    }

    listBox.innerHTML = products.map((product, index) => `
      <div class="card bg-base-100 shadow-lg border border-base-300">
        <div class="card-body">
          <div class="flex flex-col gap-3 lg:flex-row lg:items-start lg:justify-between">
            <div>
              <h2 class="card-title">Candidate ${index + 1}</h2>
              <div class="text-sm opacity-75 mt-1">
                ${escapeHtml(product.record_name)} (${escapeHtml(product.record_id)}) |
                ${product.product_start}-${product.product_end} |
                ${product.product_length} bp
              </div>
            </div>
            <form method="POST" action="${escapeHtml(downloadUrl)}">
              <input type="hidden" name="csrfmiddlewaretoken" value="${escapeHtml(getCookie("csrftoken"))}">
              <input type="hidden" name="product_sequence" value="${escapeHtml(product.product_sequence)}">
              <input type="hidden" name="pair_index" value="${escapeHtml(`${primerPair.id || "pair"}_${index + 1}`)}">
              <button type="submit" class="btn btn-outline btn-sm">Download FASTA</button>
            </form>
          </div>

          <div class="grid grid-cols-1 md:grid-cols-3 gap-4 mt-2">
            <div class="rounded-lg border border-base-300 p-3 bg-base-200/60">
              <div class="text-xs uppercase tracking-wide opacity-60">Forward binding</div>
              <div class="font-semibold mt-1">${product.forward_start}-${product.forward_end}</div>
              <div class="text-sm opacity-75 mt-1">Mismatches: ${product.forward_mismatches}</div>
            </div>
            <div class="rounded-lg border border-base-300 p-3 bg-base-200/60">
              <div class="text-xs uppercase tracking-wide opacity-60">Reverse binding</div>
              <div class="font-semibold mt-1">${product.reverse_start}-${product.reverse_end}</div>
              <div class="text-sm opacity-75 mt-1">Mismatches: ${product.reverse_mismatches}</div>
            </div>
            <div class="rounded-lg border border-base-300 p-3 bg-base-200/60">
              <div class="text-xs uppercase tracking-wide opacity-60">Primer pair</div>
              <div class="font-semibold mt-1">${escapeHtml(primerPair.name || "")}</div>
              <div class="text-sm opacity-75 mt-1">
                ${escapeHtml(primerPair.forward_name || "")} / ${escapeHtml(primerPair.reverse_name || "")}
              </div>
            </div>
          </div>

          <div class="rounded-lg border border-base-300 bg-base-200 p-4 mt-2">
            <div class="text-xs uppercase tracking-wide opacity-60">Product sequence</div>
            <pre class="mt-2 whitespace-pre-wrap break-all font-mono text-sm">${escapeHtml(product.product_sequence)}</pre>
          </div>
        </div>
      </div>
    `).join("");
  }

  async function pollStatus(taskId) {
    const response = await fetch(
      statusUrlTemplate.replace("TASK_ID", taskId),
      { headers: { Accept: "application/json" } }
    );
    if (!response.ok) {
      throw new Error("Unable to retrieve task status.");
    }
    const data = await response.json();

    if (data.state === "SUCCESS") {
      statusBox.classList.add("hidden");
      renderProducts(data.result || {});
      return;
    }

    if (data.state === "FAILURE") {
      statusBox.classList.add("hidden");
      showError(data.error || "Task failed.");
      return;
    }

    window.setTimeout(() => {
      void pollStatus(taskId);
    }, 2000);
  }

  form.addEventListener("submit", async (event) => {
    event.preventDefault();
    resetUI();
    statusBox.classList.remove("hidden");

    try {
      const response = await fetch(asyncUrl, {
        method: "POST",
        headers: { "X-CSRFToken": getCookie("csrftoken") },
        body: new FormData(form),
      });
      const data = await response.json().catch(() => ({}));

      if (!response.ok) {
        statusBox.classList.add("hidden");
        showError(data.error || "Unable to start async analysis.");
        return;
      }

      if (!data.task_id) {
        statusBox.classList.add("hidden");
        showError("Task ID missing from response.");
        return;
      }

      void pollStatus(data.task_id);
    } catch (error) {
      statusBox.classList.add("hidden");
      showError(error.message || "Unexpected error.");
    }
  });
});
