document.addEventListener("DOMContentLoaded", () => {
  const rowsNode = document.getElementById("primer-rows");
  const columnsNode = document.getElementById("primer-columns");
  const previewBody = document.getElementById("primer-preview-body");
  const editedRowsInput = document.getElementById("edited-rows");

  if (!rowsNode || !columnsNode || !previewBody || !editedRowsInput) {
    return;
  }

  const rows = JSON.parse(rowsNode.textContent || "[]");
  const columns = JSON.parse(columnsNode.textContent || "[]");
  const nameSelect = document.getElementById("id_name_column");
  const sequenceSelect = document.getElementById("id_sequence_column");
  const overhangSelect = document.getElementById("id_overhang_column");
  const submitButton = document.querySelector("form button[name='map_columns']");

  if (!columns.length || !nameSelect || !sequenceSelect || !overhangSelect) {
    return;
  }

  const buildCell = (content, className = "") => {
    const cell = document.createElement("td");
    if (className) {
      cell.className = className;
    }
    cell.appendChild(content);
    return cell;
  };

  const storePreviewRows = () => {
    const payload = Array.from(previewBody.querySelectorAll("tr")).map((row) => {
      const name = row.querySelector("td:nth-child(1) input")?.value ?? "";
      const sequence = row.querySelector("td:nth-child(2) input")?.value ?? "";
      const overhang = row.querySelector("td:nth-child(3) input")?.value ?? "";
      return { name, sequence, overhang };
    });
    editedRowsInput.value = JSON.stringify(payload);
  };

  const renderPreview = () => {
    previewBody.innerHTML = "";
    const nameKey = nameSelect.value;
    const sequenceKey = sequenceSelect.value;
    const overhangKey = overhangSelect.value;

    rows.forEach((row, index) => {
      const tr = document.createElement("tr");
      tr.dataset.index = String(index);

      const nameInput = document.createElement("input");
      nameInput.type = "text";
      nameInput.className = "input input-bordered input-sm w-full";
      nameInput.value = row[nameKey] ?? "";

      const sequenceText = document.createElement("span");
      sequenceText.className = "font-mono break-all";
      sequenceText.textContent = row[sequenceKey] ?? "";

      const sequenceInput = document.createElement("input");
      sequenceInput.type = "hidden";
      sequenceInput.value = row[sequenceKey] ?? "";

      const overhangInput = document.createElement("input");
      overhangInput.type = "text";
      overhangInput.className = "input input-bordered input-sm w-full";
      overhangInput.placeholder = "Optional";
      overhangInput.value = overhangKey ? (row[overhangKey] ?? "") : "";

      const removeButton = document.createElement("button");
      removeButton.type = "button";
      removeButton.className = "btn btn-ghost btn-sm text-error";
      removeButton.textContent = "Remove";
      removeButton.addEventListener("click", () => {
        tr.remove();
        storePreviewRows();
      });

      const sequenceWrapper = document.createElement("div");
      sequenceWrapper.appendChild(sequenceText);
      sequenceWrapper.appendChild(sequenceInput);

      tr.appendChild(buildCell(nameInput));
      tr.appendChild(buildCell(sequenceWrapper));
      tr.appendChild(buildCell(overhangInput));
      tr.appendChild(buildCell(removeButton, "text-center"));
      previewBody.appendChild(tr);
    });

    storePreviewRows();
  };

  nameSelect.addEventListener("change", renderPreview);
  sequenceSelect.addEventListener("change", renderPreview);
  overhangSelect.addEventListener("change", renderPreview);
  submitButton?.addEventListener("click", storePreviewRows);

  renderPreview();
});
