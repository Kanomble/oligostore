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
  const saveUrl = config.dataset.saveUrl;
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

    const exportFields = (product, index) => `
      <input type="hidden" name="product_sequence" value="${escapeHtml(product.product_sequence)}">
      <input type="hidden" name="pair_index" value="${escapeHtml(`${primerPair.id || "pair"}_${index + 1}`)}">
      <input type="hidden" name="sequence_file_id" value="${escapeHtml(sequenceFile.id || "")}">
      <input type="hidden" name="record_id" value="${escapeHtml(product.record_id || "")}">
      <input type="hidden" name="product_start" value="${escapeHtml(product.product_start)}">
      <input type="hidden" name="product_end" value="${escapeHtml(product.product_end)}">
      <input type="hidden" name="forward_overhang_sequence" value="${escapeHtml(product.forward_overhang_sequence || "")}">
      <input type="hidden" name="reverse_overhang_sequence" value="${escapeHtml(product.reverse_overhang_sequence || "")}">
      <input type="hidden" name="wraps_origin" value="${escapeHtml(product.wraps_origin ? "true" : "false")}">
    `;

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
              ${(product.is_circular_record || product.wraps_origin) ? `
                <div class="text-xs opacity-70 mt-1">
                  ${product.wraps_origin ? "Circular amplicon crossing the record origin." : "Detected on a circular record."}
                </div>
              ` : ""}
            </div>
            <div class="flex flex-wrap gap-2">
              <form method="POST" action="${escapeHtml(downloadUrl)}">
                <input type="hidden" name="csrfmiddlewaretoken" value="${escapeHtml(getCookie("csrftoken"))}">
                <input type="hidden" name="export_format" value="fasta">
                ${exportFields(product, index)}
                <button type="submit" class="btn btn-outline btn-sm">Download FASTA</button>
              </form>
              <form method="POST" action="${escapeHtml(downloadUrl)}">
                <input type="hidden" name="csrfmiddlewaretoken" value="${escapeHtml(getCookie("csrftoken"))}">
                <input type="hidden" name="export_format" value="genbank">
                ${exportFields(product, index)}
                <button type="submit" class="btn btn-outline btn-sm">Download GenBank</button>
              </form>
            </div>
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
            <div class="flex items-center justify-between gap-3">
              <div class="text-xs uppercase tracking-wide opacity-60">Product sequence</div>
              <div class="text-xs opacity-60">${product.product_sequence.length} nt</div>
            </div>
            <textarea class="textarea textarea-bordered mt-2 w-full font-mono text-sm" rows="6" readonly>${escapeHtml(product.product_sequence)}</textarea>
            <div class="mt-3 flex flex-col gap-2 sm:flex-row sm:items-end">
              <label class="form-control flex-1">
                <span class="label py-1"><span class="label-text text-xs">Save name</span></span>
                <input
                  type="text"
                  class="input input-bordered input-sm"
                  data-save-product-name="${index}"
                  value="${escapeHtml(`${primerPair.name || "PCR_Product"}_${product.record_id}_${product.product_start}-${product.product_end}`)}"
                />
              </label>
              <button
                type="button"
                class="btn btn-primary btn-sm"
                data-save-product-index="${index}"
              >
                Save PCR product
              </button>
            </div>
            <div class="text-xs opacity-70" data-save-product-status="${index}"></div>
          </div>
        </div>
      </div>
    `).join("");
  }

  async function saveProduct(productIndex) {
    const resultProducts = window.__primerPairProductsResult?.products || [];
    const primerPair = window.__primerPairProductsResult?.primer_pair || {};
    const sequenceFile = window.__primerPairProductsResult?.sequence_file || {};
    const product = resultProducts[productIndex];
    const nameInput = listBox.querySelector(`[data-save-product-name="${productIndex}"]`);
    const statusEl = listBox.querySelector(`[data-save-product-status="${productIndex}"]`);
    const button = listBox.querySelector(`[data-save-product-index="${productIndex}"]`);

    if (!product || !nameInput || !statusEl || !button) {
      return;
    }

    button.disabled = true;
    statusEl.textContent = "Saving PCR product...";
    statusEl.className = "text-xs opacity-70";

    try {
      const response = await fetch(saveUrl, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          "X-CSRFToken": getCookie("csrftoken"),
        },
        body: JSON.stringify({
          primer_pair_id: primerPair.id,
          sequence_file_id: sequenceFile.id,
          name: nameInput.value.trim(),
          record_id: product.record_id,
          product_start: product.product_start,
          product_end: product.product_end,
          product_sequence: product.product_sequence,
          wraps_origin: product.wraps_origin,
        }),
      });
      const payload = await response.json().catch(() => ({}));
      if (!response.ok) {
        throw new Error(payload.error || "Unable to save PCR product.");
      }
      nameInput.value = payload.pcr_product?.name || nameInput.value;
      statusEl.textContent = `Saved as ${payload.pcr_product?.name || "PCR product"}.`;
      statusEl.className = "text-xs text-success";
    } catch (error) {
      statusEl.textContent = error.message || "Unable to save PCR product.";
      statusEl.className = "text-xs text-error";
    } finally {
      button.disabled = false;
    }
  }

  async function pollStatus(jobId) {
    const response = await fetch(
      statusUrlTemplate.replace("TASK_ID", jobId),
      { headers: { Accept: "application/json" } }
    );
    if (!response.ok) {
      throw new Error("Unable to retrieve task status.");
    }
    const data = await response.json();

    if (data.state === "SUCCESS") {
      statusBox.classList.add("hidden");
      window.__primerPairProductsResult = data.result || {};
      renderProducts(data.result || {});
      return;
    }

    if (data.state === "FAILURE") {
      statusBox.classList.add("hidden");
      showError(data.error || "Task failed.");
      return;
    }

    window.setTimeout(() => {
      void pollStatus(jobId);
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

      const jobId = data.job_id || data.task_id;
      if (!jobId) {
        statusBox.classList.add("hidden");
        showError("Job ID missing from response.");
        return;
      }

      void pollStatus(jobId);
    } catch (error) {
      statusBox.classList.add("hidden");
      showError(error.message || "Unexpected error.");
    }
  });

  listBox.addEventListener("click", (event) => {
    const button = event.target.closest("[data-save-product-index]");
    if (!button) {
      return;
    }
    void saveProduct(Number(button.dataset.saveProductIndex));
  });
});
