(function () {
  const fragmentAnchors = new Map();

  function getCloningPreviewForm() {
    return document.getElementById("cloning-assembly-preview-form");
  }

  function setCloningVisualFeedback(message, tone) {
    const feedback = document.querySelector("[data-cloning-visual-feedback]");
    if (!feedback) {
      return;
    }
    feedback.textContent = message || "";
    feedback.classList.remove("text-error", "text-success", "text-warning");
    if (tone) {
      feedback.classList.add(tone);
    }
  }

  function ensureSelectOption(select, value) {
    if (!select || !value) {
      return false;
    }
    for (const option of select.options) {
      if (option.value === value) {
        return true;
      }
    }
    return false;
  }

  function setSelectValue(select, value) {
    if (!select) {
      return false;
    }
    if (value && !ensureSelectOption(select, value)) {
      return false;
    }
    select.value = value || "";
    select.dispatchEvent(new Event("change", { bubbles: true }));
    return true;
  }

  function setFormLoading(form, message) {
    if (!form) {
      return;
    }
    form.classList.add("is-loading");
    form.setAttribute("aria-busy", "true");
    form.querySelectorAll('button[type="submit"]').forEach((button) => {
      button.classList.add("loading");
      button.disabled = true;
    });
    if (message) {
      setCloningVisualFeedback(message, "text-warning");
    }
  }

  function submitForm(form, message) {
    if (!form) {
      return;
    }
    setFormLoading(form, message || "Updating assembly preview...");
    if (typeof form.requestSubmit === "function") {
      form.requestSubmit();
      return;
    }
    HTMLFormElement.prototype.submit.call(form);
  }

  function refreshVisualEnzymeSelection(form) {
    if (!form) {
      return;
    }
    const left = form.elements.left_enzyme ? form.elements.left_enzyme.value : "";
    const right = form.elements.right_enzyme ? form.elements.right_enzyme.value : "";
    document.querySelectorAll("[data-cloning-enzyme-picker]").forEach((site) => {
      const enzymeName = site.dataset.enzymeName || "";
      site.classList.toggle("is-selected", Boolean(enzymeName && (enzymeName === left || enzymeName === right)));
    });
  }

  function submitPreviewWhenReady(form) {
    if (!form) {
      return;
    }
    const left = form.elements.left_enzyme ? form.elements.left_enzyme.value : "";
    const right = form.elements.right_enzyme ? form.elements.right_enzyme.value : "";
    refreshVisualEnzymeSelection(form);
    if (left && right) {
      submitForm(form, "Updating assembly preview...");
      return;
    }
    setCloningVisualFeedback("Select another restriction site to complete the enzyme pair.", "text-warning");
  }

  function clearFragmentAnchors(role) {
    const roles = role ? [role] : Array.from(fragmentAnchors.keys());
    roles.forEach((mapRole) => {
      const previous = fragmentAnchors.get(mapRole);
      if (previous && previous.element) {
        previous.element.classList.remove("is-fragment-anchor");
      }
      fragmentAnchors.delete(mapRole);
    });
  }

  function siteData(trigger) {
    const map = trigger.closest("[data-cloning-sequence-map]");
    const position = Number.parseInt(trigger.dataset.sitePosition || "", 10);
    return {
      element: trigger,
      role: trigger.dataset.mapRole || (map ? map.dataset.mapRole : ""),
      mapShape: map ? map.dataset.mapShape || "linear" : "linear",
      enzymeName: trigger.dataset.enzymeName || "",
      siteId: trigger.dataset.siteId || "",
      position: Number.isFinite(position) ? position : null,
    };
  }

  function fragmentData(trigger) {
    const start = Number.parseInt(trigger.dataset.fragmentStart || "", 10);
    const end = Number.parseInt(trigger.dataset.fragmentEnd || "", 10);
    return {
      element: trigger,
      role: trigger.dataset.mapRole || "",
      fieldName: trigger.dataset.fragmentField || "",
      index: trigger.dataset.fragmentIndex || "",
      start: Number.isFinite(start) ? start : null,
      end: Number.isFinite(end) ? end : null,
      wrapsOrigin: trigger.dataset.fragmentWrapsOrigin === "1",
    };
  }

  function selectedOverlayEnzymes(form) {
    const select = form && form.elements.selected_enzymes;
    if (!select || !select.options) {
      return [];
    }
    return Array.from(select.options)
      .filter((option) => option.selected && option.value)
      .map((option) => option.value);
  }

  function digestRowsCanRepresentSameEnzymeFragments(form, enzymeName) {
    const overlays = selectedOverlayEnzymes(form);
    if (overlays.length === 0) {
      return true;
    }
    return overlays.length === 1 && overlays[0] === enzymeName;
  }

  function fragmentMatchesPair(fragment, firstSite, secondSite) {
    if (fragment.start === null || fragment.end === null || firstSite.position === null || secondSite.position === null) {
      return false;
    }
    if (firstSite.mapShape === "circular") {
      return (
        fragment.start === firstSite.position &&
        fragment.end === secondSite.position &&
        fragment.wrapsOrigin === (secondSite.position <= firstSite.position)
      );
    }
    return (
      fragment.start === Math.min(firstSite.position, secondSite.position) &&
      fragment.end === Math.max(firstSite.position, secondSite.position)
    );
  }

  function matchingFragmentForPair(form, firstSite, secondSite) {
    const regionCandidates = Array.from(
      document.querySelectorAll("[data-cloning-region-picker]")
    )
      .filter((candidate) => candidate.dataset.mapRole === firstSite.role)
      .map(fragmentData);
    const regionMatch = regionCandidates.find((fragment) => fragmentMatchesPair(fragment, firstSite, secondSite));
    if (regionMatch) {
      return regionMatch;
    }

    if (!digestRowsCanRepresentSameEnzymeFragments(form, firstSite.enzymeName)) {
      return null;
    }
    return Array.from(
      document.querySelectorAll("[data-cloning-digest-fragment]")
    )
      .filter((candidate) => candidate.dataset.mapRole === firstSite.role)
      .map(fragmentData)
      .find((fragment) => fragmentMatchesPair(fragment, firstSite, secondSite)) || null;
  }

  function handleShiftClickSite(trigger, form) {
    const current = siteData(trigger);
    if (!current.role || !current.enzymeName || current.position === null) {
      setCloningVisualFeedback("That restriction site is missing fragment-selection metadata.", "text-error");
      return;
    }

    const previous = fragmentAnchors.get(current.role);
    if (!previous || previous.siteId === current.siteId) {
      clearFragmentAnchors(current.role);
      fragmentAnchors.set(current.role, current);
      trigger.classList.add("is-fragment-anchor");
      setCloningVisualFeedback(`Fragment anchor set on ${current.role} ${current.enzymeName}. Shift-click the second site on the same map.`, "text-warning");
      return;
    }

    if (previous.enzymeName !== current.enzymeName) {
      clearFragmentAnchors(current.role);
      fragmentAnchors.set(current.role, current);
      trigger.classList.add("is-fragment-anchor");
      setCloningVisualFeedback("Fragment selection currently maps to same-enzyme digest fragments. Anchor reset to this site.", "text-warning");
      return;
    }

    const fragment = matchingFragmentForPair(form, previous, current);
    if (!fragment || !fragment.fieldName || !fragment.index) {
      clearFragmentAnchors(current.role);
      setCloningVisualFeedback("No selectable digest fragment matches those two sites. Use adjacent same-enzyme cut sites or select the fragment from the map.", "text-error");
      return;
    }

    const field = form.elements[fragment.fieldName];
    const left = form.elements.left_enzyme;
    const right = form.elements.right_enzyme;
    if (!field || !left || !right || !ensureSelectOption(left, current.enzymeName) || !ensureSelectOption(right, current.enzymeName)) {
      clearFragmentAnchors(current.role);
      setCloningVisualFeedback("That fragment cannot be submitted with the current enzyme choices.", "text-error");
      return;
    }

    setSelectValue(left, current.enzymeName);
    setSelectValue(right, current.enzymeName);
    field.value = fragment.index;
    clearFragmentAnchors(current.role);
    submitForm(form, `Selecting ${current.role} fragment ${fragment.index}...`);
  }

  document.addEventListener("click", function (event) {
    const regionTrigger = event.target.closest("[data-cloning-region-picker], [data-cloning-fragment-picker]");
    if (regionTrigger) {
      const form = getCloningPreviewForm();
      if (!form) {
        return;
      }
      const field = form.elements[regionTrigger.dataset.fragmentField];
      if (!field) {
        return;
      }
      const nextValue = regionTrigger.dataset.fragmentIndex || "";
      field.value = field.value === nextValue ? "" : nextValue;
      submitForm(form, "Updating selected fragment...");
      return;
    }

    const enzymeTrigger = event.target.closest("[data-cloning-enzyme-picker]");
    if (!enzymeTrigger) {
      return;
    }

    const form = getCloningPreviewForm();
    if (!form) {
      return;
    }

    if (event.shiftKey) {
      handleShiftClickSite(enzymeTrigger, form);
      return;
    }

    clearFragmentAnchors();
    const enzymeName = enzymeTrigger.dataset.enzymeName || "";
    const left = form.elements.left_enzyme;
    const right = form.elements.right_enzyme;
    if (!enzymeName || !left || !right || !ensureSelectOption(left, enzymeName) || !ensureSelectOption(right, enzymeName)) {
      setCloningVisualFeedback("That restriction site is not available for the selected vector.", "text-error");
      return;
    }

    if (left.value === enzymeName && right.value === enzymeName) {
      setSelectValue(right, "");
      setCloningVisualFeedback("Right enzyme cleared. Click another site for directional cloning or click Preview Assembly after choosing a second enzyme.", "text-warning");
      refreshVisualEnzymeSelection(form);
      return;
    }
    if (!left.value) {
      setSelectValue(left, enzymeName);
      if (!right.value) {
        setSelectValue(right, enzymeName);
      }
      submitPreviewWhenReady(form);
      return;
    }
    if (left.value === enzymeName && right.value !== enzymeName) {
      setSelectValue(right, enzymeName);
      submitPreviewWhenReady(form);
      return;
    }
    if (right.value === enzymeName && left.value !== enzymeName) {
      setSelectValue(left, enzymeName);
      submitPreviewWhenReady(form);
      return;
    }
    if (left.value === right.value) {
      setSelectValue(right, enzymeName);
      submitPreviewWhenReady(form);
      return;
    }
    setSelectValue(right, enzymeName);
    submitPreviewWhenReady(form);
  });

  document.addEventListener("keydown", function (event) {
    if (event.key !== "Enter" && event.key !== " ") {
      return;
    }
    const trigger = event.target.closest("[data-cloning-enzyme-picker], [data-cloning-region-picker]");
    if (!trigger) {
      return;
    }
    event.preventDefault();
    trigger.dispatchEvent(new MouseEvent("click", { bubbles: true, shiftKey: event.shiftKey }));
  });

  document.addEventListener("submit", function (event) {
    if (event.target && event.target.id === "cloning-assembly-preview-form") {
      setFormLoading(event.target, "Updating assembly preview...");
    }
  });

  refreshVisualEnzymeSelection(getCloningPreviewForm());
})();
