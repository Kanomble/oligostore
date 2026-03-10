document.addEventListener("DOMContentLoaded", () => {
  const selectionContainer = document.getElementById("primerpair-selection");
  if (!selectionContainer) {
    return;
  }

  const forwardInput = document.getElementById("id_forward_primer");
  const reverseInput = document.getElementById("id_reverse_primer");
  const forwardSummary = document.getElementById("selected-forward");
  const reverseSummary = document.getElementById("selected-reverse");
  const stepLabel = document.getElementById("selection-step");
  const createButton = document.getElementById("create-primer-pair");
  const actionButtons = selectionContainer.querySelectorAll("[data-action]");
  const primerButtons = document.querySelectorAll("[data-primer-select]");
  const primerRows = document.querySelectorAll("[data-primer-row]");
  const storageKey = "primerpair-create-selection";

  let forwardId = selectionContainer.dataset.selectedForwardId || "";
  let forwardName = selectionContainer.dataset.selectedForwardName || "";
  let reverseId = selectionContainer.dataset.selectedReverseId || "";
  let reverseName = selectionContainer.dataset.selectedReverseName || "";
  let activeStep = "forward";

  function loadStoredSelection() {
    try {
      const stored = sessionStorage.getItem(storageKey);
      return stored ? JSON.parse(stored) : null;
    } catch {
      return null;
    }
  }

  function saveStoredSelection() {
    if (!forwardId && !reverseId) {
      clearStoredSelection();
      return;
    }
    try {
      sessionStorage.setItem(storageKey, JSON.stringify({
        forwardId,
        forwardName,
        reverseId,
        reverseName,
      }));
    } catch {
      return;
    }
  }

  function clearStoredSelection() {
    try {
      sessionStorage.removeItem(storageKey);
    } catch {
      return;
    }
  }

  if (!forwardId && !reverseId) {
    const storedSelection = loadStoredSelection();
    if (storedSelection) {
      forwardId = storedSelection.forwardId || "";
      forwardName = storedSelection.forwardName || "";
      reverseId = storedSelection.reverseId || "";
      reverseName = storedSelection.reverseName || "";
    }
  }

  if (forwardId) {
    activeStep = "reverse";
  }

  function setActiveStep(step) {
    activeStep = step;
    updateUI();
  }

  function updateStepLabel() {
    if (!forwardId) {
      stepLabel.textContent = "Select a forward primer.";
      return;
    }
    if (!reverseId) {
      stepLabel.textContent = "Select a reverse primer.";
      return;
    }
    stepLabel.textContent = activeStep === "forward"
      ? "Change forward primer if needed."
      : "Change reverse primer if needed.";
  }

  function updateSummary() {
    forwardSummary.textContent = forwardName ? `${forwardName} (ID ${forwardId})` : "None selected";
    reverseSummary.textContent = reverseName ? `${reverseName} (ID ${reverseId})` : "None selected";
  }

  function updateHiddenInputs() {
    forwardInput.value = forwardId || "";
    reverseInput.value = reverseId || "";
  }

  function updateCreateButton() {
    createButton.disabled = !(forwardId && reverseId);
  }

  function updateRowStyles() {
    primerRows.forEach((row) => {
      const primerId = row.dataset.primerId;
      row.classList.toggle("bg-primary/10", primerId === forwardId);
      row.classList.toggle("bg-secondary/10", primerId === reverseId);
    });
  }

  function updateSelectButtons() {
    primerButtons.forEach((button) => {
      const primerId = button.dataset.primerId;
      button.disabled = activeStep === "reverse" && primerId === forwardId;
      if (primerId === forwardId) {
        button.textContent = "Forward";
        return;
      }
      if (primerId === reverseId) {
        button.textContent = "Reverse";
        return;
      }
      button.textContent = activeStep === "forward" ? "Select Forward" : "Select Reverse";
    });
  }

  function updateUI() {
    updateStepLabel();
    updateSummary();
    updateHiddenInputs();
    updateCreateButton();
    updateRowStyles();
    updateSelectButtons();
    saveStoredSelection();
  }

  function clearForward() {
    forwardId = "";
    forwardName = "";
    reverseId = "";
    reverseName = "";
    activeStep = "forward";
    clearStoredSelection();
    updateUI();
  }

  function clearReverse() {
    reverseId = "";
    reverseName = "";
    activeStep = "reverse";
    saveStoredSelection();
    updateUI();
  }

  primerButtons.forEach((button) => {
    button.addEventListener("click", () => {
      const primerId = button.dataset.primerId;
      const primerName = button.dataset.primerName;

      if (activeStep === "forward") {
        forwardId = primerId;
        forwardName = primerName;
        if (reverseId === primerId) {
          reverseId = "";
          reverseName = "";
        }
        activeStep = "reverse";
        updateUI();
        return;
      }

      if (primerId === forwardId) {
        return;
      }
      reverseId = primerId;
      reverseName = primerName;
      updateUI();
    });
  });

  actionButtons.forEach((button) => {
    button.addEventListener("click", () => {
      const action = button.dataset.action;
      if (action === "change-forward") {
        setActiveStep("forward");
        return;
      }
      if (action === "change-reverse") {
        setActiveStep("reverse");
        return;
      }
      if (action === "clear-forward") {
        clearForward();
        return;
      }
      if (action === "clear-reverse") {
        clearReverse();
      }
    });
  });

  updateUI();
});
