document.addEventListener("DOMContentLoaded", () => {
  const form = document.getElementById("primer-binding-form");
  const config = document.getElementById("primer-binding-config");
  const statusBox = document.getElementById("primer-binding-status");
  const errorBox = document.getElementById("primer-binding-error");
  const resultsBox = document.getElementById("primer-binding-results");

  if (!form || !config || !statusBox || !errorBox || !resultsBox) {
    return;
  }

  const errorMessage = errorBox.querySelector(".alert");
  const emptyBox = document.getElementById("primer-binding-empty");
  const tableBox = document.getElementById("primer-binding-table");
  const tableBody = tableBox ? tableBox.querySelector("tbody") : null;
  const asyncUrl = config.dataset.asyncUrl;
  const statusUrlTemplate = config.dataset.statusUrlTemplate;

  function getCookie(name) {
    const value = `; ${document.cookie}`;
    const parts = value.split(`; ${name}=`);
    if (parts.length === 2) {
      return parts.pop().split(";").shift();
    }
    return "";
  }

  function resetUI() {
    statusBox.classList.add("hidden");
    errorBox.classList.add("hidden");
    resultsBox.classList.add("hidden");
    emptyBox.classList.add("hidden");
    tableBox.classList.add("hidden");
    if (tableBody) {
      tableBody.innerHTML = "";
    }
  }

  function showError(message) {
    errorMessage.textContent = message;
    errorBox.classList.remove("hidden");
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
      resultsBox.classList.remove("hidden");
      const hits = data.result || [];
      if (!hits.length) {
        emptyBox.classList.remove("hidden");
        return;
      }

      tableBox.classList.remove("hidden");
      hits.forEach((hit) => {
        const row = document.createElement("tr");
        row.innerHTML = `
          <td>${hit.record_id}</td>
          <td>${hit.start}</td>
          <td>${hit.end}</td>
          <td>${hit.strand}</td>
          <td>${hit.mismatches}</td>
        `;
        tableBody.appendChild(row);
      });
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

    const formData = new FormData(form);
    statusBox.classList.remove("hidden");

    try {
      const response = await fetch(asyncUrl, {
        method: "POST",
        headers: { "X-CSRFToken": getCookie("csrftoken") },
        body: formData,
      });

      if (!response.ok) {
        statusBox.classList.add("hidden");
        showError("Unable to start async analysis.");
        return;
      }

      const data = await response.json();
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
