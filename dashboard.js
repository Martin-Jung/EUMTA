const YEARS = ["<2000", "2000-2006", "2006-2012", "2012-2018", "2018-2024"];
const TARGET_LABELS = { loglinear: "Log-linear", flat: "Flat", extinction: "Extinction" };
const COLORS = {
  combined: "var(--chart-combined)", natura2000: "var(--chart-natura2000)", cdda: "var(--chart-cdda)",
  art12: "var(--chart-art12)", art17: "var(--chart-art17)", Habitats: "var(--chart-habitats)", Species: "var(--chart-species)"
};

/** Escapes dynamic text before it is inserted into chart markup.
 * @param {*} value Value to escape.
 * @returns {string} HTML-safe text.
 */
function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, (character) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", "\"": "&quot;", "'": "&#39;"
  })[character]);
}

/** Builds the legend markup for a set of chart series.
 * @param {Array<{label: string, color: string}>} series Legend labels and colors.
 * @returns {string} Legend HTML.
 */
function renderLegend(series) {
  return `<div class="chart-legend">${series.map((item) => `<span class="legend-item"><span class="legend-swatch" style="background:${item.color}"></span>${escapeHtml(item.label)}</span>`).join("")}</div>`;
}

/**
 * Renders an SVG line chart and its legend into a container.
 * @param {HTMLElement} container Chart host element.
 * @param {Object[]} rows Source rows consumed by the series.
 * @param {Array<{label: string, color: string, filter: Function, value: Function}>} series Lines to render.
 * @param {{width?: number, height?: number, min?: number, max?: number, ticks?: number[], threshold?: number, formatTick?: Function, yAxisTitle?: string, label?: string}} [options] Dimensions, scale, ticks, threshold, and accessible label.
 * @returns {void} Replaces the container contents with the chart and legend.
 */
function renderLineChart(container, rows, series, options = {}) {
  const width = options.width || 900;
  const height = options.height || 320;
  const margin = { top: 20, right: 78, bottom: 58, left: 72 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const min = options.min ?? 0;
  const max = options.max ?? 1;
  const ticks = options.ticks || [0, 0.25, 0.5, 0.75, 1];
  const x = (index) => margin.left + (index / (YEARS.length - 1)) * plotWidth;
  const y = (value) => margin.top + ((max - value) / (max - min)) * plotHeight;
  const formatTick = options.formatTick || ((value) => Number(value).toFixed(2));
  const parts = [`<svg class="chart-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="${escapeHtml(options.label || "MTA trend chart")}">`];
  ticks.forEach((tick) => {
    const position = y(tick);
    parts.push(`<line class="grid-line" x1="${margin.left}" y1="${position}" x2="${width - margin.right}" y2="${position}"></line>`);
    parts.push(`<text x="${margin.left - 10}" y="${position + 4}" text-anchor="end">${escapeHtml(formatTick(tick))}</text>`);
  });
  if (options.threshold !== undefined) {
    const position = y(options.threshold);
    parts.push(`<line class="target-line" x1="${margin.left}" y1="${position}" x2="${width - margin.right}" y2="${position}"></line>`);
  }
  parts.push(`<line class="axis-line" x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top + plotHeight}"></line>`);
  parts.push(`<line class="axis-line" x1="${margin.left}" y1="${margin.top + plotHeight}" x2="${width - margin.right}" y2="${margin.top + plotHeight}"></line>`);
  parts.push(`<text class="axis-title" transform="translate(17 ${margin.top + plotHeight / 2}) rotate(-90)" text-anchor="middle">${escapeHtml(options.yAxisTitle || "MTA")}</text>`);
  YEARS.forEach((year, index) => parts.push(`<text x="${x(index)}" y="${height - 27}" text-anchor="middle">${escapeHtml(year)}</text>`));
  series.forEach((item) => {
    const itemRows = rows.filter(item.filter);
    const points = YEARS.map((year, index) => {
      const row = itemRows.find((candidate) => candidate.year === year);
      return row && Number.isFinite(Number(item.value(row))) ? { x: x(index), y: y(Number(item.value(row))), year, value: Number(item.value(row)) } : null;
    }).filter(Boolean);
    if (!points.length) return;
    const path = points.map((point, index) => `${index ? "L" : "M"}${point.x},${point.y}`).join(" ");
    parts.push(`<path class="series-line" d="${path}" stroke="${item.color}"></path>`);
    points.forEach((point) => {
      const title = `${item.label} · ${point.year}: ${formatTick(point.value)}`;
      parts.push(`<circle class="series-point" cx="${point.x}" cy="${point.y}" r="4.5" fill="${item.color}"><title>${escapeHtml(title)}</title></circle>`);
    });
  });
  parts.push("</svg>");
  container.innerHTML = `${parts.join("")}${renderLegend(series)}`;
}

/** Displays the empty-state message for a target method with no data.
 * @param {HTMLElement} container Chart host element.
 * @param {string} targetOption Target method used for the message.
 * @returns {void} Replaces the container contents with the empty state.
 */
function renderTargetEmpty(container, targetOption) {
  container.innerHTML = `<div class="chart-loading chart-empty" role="status">No data available for ${escapeHtml(TARGET_LABELS[targetOption])} targets yet.</div>`;
}

/** Renders a target-dependent line chart or its empty state.
 * @param {HTMLElement} container Chart host element.
 * @param {Object[]} rows Source rows.
 * @param {Array<{label: string, color: string, filter: Function, value: Function}>} series Candidate chart lines.
 * @param {Object} options Line-chart options passed to renderLineChart.
 * @param {string} targetOption Selected target method.
 * @returns {void} Replaces the container contents with a chart or empty state.
 */
function renderTargetLineChart(container, rows, series, options, targetOption) {
  if (!rows.some((row) => series.some((item) => item.filter(row)))) {
    renderTargetEmpty(container, targetOption);
    return;
  }
  renderLineChart(container, rows, series, options);
}

/** Renders grouped country bars for feature counts or distribution coverage.
 * @param {HTMLElement} container Chart host element.
 * @param {Object[]} rows Country summary rows.
 * @param {string} metric Row field to chart: feature count or coverage proportion.
 * @returns {void} Replaces the container contents with the SVG chart and legend.
 */
function renderCountryBars(container, rows, metric) {
  const selected = rows.filter((row) => row.ctype === "combined" && Number.isFinite(Number(row[metric])));
  const countries = [...new Set(selected.map((row) => row.country))].sort((first, second) => first.localeCompare(second));
  const width = 900;
  const rowHeight = 27;
  const height = countries.length * rowHeight + 94;
  const margin = { top: 24, right: 42, bottom: 68, left: 72 };
  const plotWidth = width - margin.left - margin.right;
  const xAxisTitle = metric === "prop" ? "Share of distribution" : "Count";
  const maxValue = metric === "prop" ? 100 : Math.ceil(Math.max(...selected.map((row) => row[metric])) / 100) * 100;
  const ticks = metric === "prop" ? [0, 25, 50, 75, 100] : Array.from({ length: 5 }, (_, index) => (maxValue * index) / 4);
  const groups = ["Habitats", "Species"];
  const parts = [`<svg class="chart-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="${metric === "prop" ? "Mean distribution conserved by country" : "Feature counts by country"}">`];
  ticks.forEach((tick) => {
    const position = margin.left + (tick / maxValue) * plotWidth;
    parts.push(`<line class="grid-line" x1="${position}" y1="${margin.top - 10}" x2="${position}" y2="${height - margin.bottom + 3}"></line>`);
    parts.push(`<text x="${position}" y="${height - 43}" text-anchor="middle">${metric === "prop" ? `${Math.round(tick)}%` : Math.round(tick)}</text>`);
  });
  countries.forEach((country, index) => {
    const top = margin.top + index * rowHeight;
    parts.push(`<text x="${margin.left - 13}" y="${top + 14}" text-anchor="end">${escapeHtml(country)}</text>`);
    groups.forEach((group, groupIndex) => {
      const item = selected.find((candidate) => candidate.country === country && candidate.category === group);
      if (!item) return;
      const rawValue = Number(item[metric]);
      const value = metric === "prop" ? rawValue * 100 : rawValue;
      const barWidth = (value / maxValue) * plotWidth;
      const barY = top + (groupIndex === 0 ? 2 : 14);
      const label = metric === "prop" ? `${value.toFixed(1)}%` : Math.round(value).toLocaleString("en-US");
      parts.push(`<rect x="${margin.left}" y="${barY}" width="${barWidth}" height="8" rx="2" fill="${COLORS[group]}"><title>${escapeHtml(country)} · ${escapeHtml(group)}: ${label}</title></rect>`);
      if (barWidth > 35) parts.push(`<text x="${margin.left + barWidth + 5}" y="${barY + 8}" font-size="9">${label}</text>`);
    });
  });
  parts.push(`<line class="axis-line" x1="${margin.left}" y1="${margin.top - 10}" x2="${margin.left}" y2="${height - margin.bottom + 3}"></line>`);
  parts.push(`<line class="axis-line" x1="${margin.left}" y1="${height - margin.bottom + 3}" x2="${width - margin.right}" y2="${height - margin.bottom + 3}"></line>`);
  parts.push(`<text class="axis-title" transform="translate(17 ${(margin.top - 10 + height - margin.bottom + 3) / 2}) rotate(-90)" text-anchor="middle">Country</text>`);
  parts.push(`<text class="axis-title x-axis-title" x="${margin.left + plotWidth / 2}" y="${height - 8}" text-anchor="middle">${xAxisTitle}</text>`);
  parts.push("</svg>");
  container.innerHTML = `${parts.join("")}${renderLegend(groups.map((label) => ({ label, color: COLORS[label] })))}`;
  const svg = container.querySelector("svg");
  const title = svg.querySelector(".x-axis-title");
  const alignAxisTitle = () => {
    const scale = width / svg.getBoundingClientRect().width;
    const visibleStart = container.scrollLeft * scale;
    const visibleEnd = visibleStart + container.clientWidth * scale;
    const axisStart = Math.max(margin.left, visibleStart);
    const axisEnd = Math.min(width - margin.right, visibleEnd);
    title.setAttribute("x", String((axisStart + axisEnd) / 2));
  };
  alignAxisTitle();
  container.onscroll = alignAxisTitle;
}

/** Bins feature target values and renders their distribution as an SVG histogram.
 * @param {HTMLElement} container Chart host element.
 * @param {Object[]} rows Target rows to bin.
 * @param {string} targetOption Selected method for an empty-state message.
 * @returns {void} Replaces the container contents with a histogram or empty state.
 */
function renderHistogram(container, rows, targetOption) {
  if (!rows.length) {
    renderTargetEmpty(container, targetOption);
    return;
  }
  const values = rows.map((row) => Number(row.target_relative)).filter((value) => Number.isFinite(value) && value >= 0 && value <= 1);
  const binCount = 30;
  const bins = Array.from({ length: binCount }, () => 0);
  values.forEach((value) => { bins[Math.min(binCount - 1, Math.floor(value * binCount))] += 1; });
  const width = 900;
  const height = 300;
  const margin = { top: 18, right: 20, bottom: 48, left: 78 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const maxCount = Math.ceil(Math.max(...bins) / 100) * 100 || 1;
  const parts = [`<svg class="chart-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="Distribution of feature targets">`];
  [0, 0.25, 0.5, 0.75, 1].forEach((fraction) => {
    const position = margin.top + plotHeight * (1 - fraction);
    parts.push(`<line class="grid-line" x1="${margin.left}" y1="${position}" x2="${width - margin.right}" y2="${position}"></line>`);
    parts.push(`<text x="${margin.left - 10}" y="${position + 4}" text-anchor="end">${Math.round(maxCount * fraction)}</text>`);
  });
  bins.forEach((count, index) => {
    const barWidth = plotWidth / binCount;
    const barHeight = (count / maxCount) * plotHeight;
    const barX = margin.left + index * barWidth + 1;
    const barY = margin.top + plotHeight - barHeight;
    parts.push(`<rect x="${barX}" y="${barY}" width="${Math.max(1, barWidth - 2)}" height="${barHeight}" fill="${COLORS.art17}"><title>${(index / binCount).toFixed(2)}-${((index + 1) / binCount).toFixed(2)}: ${count} features</title></rect>`);
  });
  [0, 0.25, 0.5, 0.75, 1].forEach((fraction) => {
    const position = margin.left + plotWidth * fraction;
    parts.push(`<text x="${position}" y="${height - 17}" text-anchor="middle">${Math.round(fraction * 100)}%</text>`);
  });
  parts.push(`<line class="axis-line" x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${margin.top + plotHeight}"></line>`);
  parts.push(`<line class="axis-line" x1="${margin.left}" y1="${margin.top + plotHeight}" x2="${width - margin.right}" y2="${margin.top + plotHeight}"></line>`);
  parts.push(`<text class="axis-title" transform="translate(17 ${margin.top + plotHeight / 2}) rotate(-90)" text-anchor="middle">Number of features</text>`);
  parts.push("</svg>");
  container.innerHTML = parts.join("");
}

/** Updates summary values and renders every dashboard chart for the selected target method.
 * @param {{overall: Object[], statistics: Object[], countries: Object[], targets: Object[], areas: Object[], art12: Object[], art17Species: Object[], art17Habitats: Object[]}} data Dashboard datasets.
 * @param {string} [targetOption="loglinear"] Target method used to filter MTA data.
 * @returns {void} Updates summary and chart elements in the document.
 */
function renderDashboard(data, targetOption = "loglinear") {
  const overall = data.overall.filter((row) => row.option === targetOption);
  const combined = overall.filter((row) => row.ctype === "combined" && row.category === "All");
  const latest = combined.find((row) => row.year === "2018-2024") || combined[combined.length - 1];
  const previous = combined.find((row) => row.year === "2012-2018");
  const featureCount = [0, 2, 4].reduce((total, index) => total + Number(data.statistics[index]?.nr_features || 0), 0);
  document.querySelector("#target-status").textContent = overall.length ? "" : `${TARGET_LABELS[targetOption]} targets are not computed yet.`;
  document.querySelector("#latest-mta").textContent = latest ? Number(latest.y).toFixed(3) : "—";
  document.querySelector("#change-mta").textContent = latest && previous ? `${latest.y >= previous.y ? "+" : ""}${(latest.y - previous.y).toFixed(3)}` : "—";
  document.querySelector("#feature-count").textContent = featureCount.toLocaleString("en-US");

  renderTargetLineChart(document.querySelector("#overall-chart"), overall, [
    { label: "Combined areas", color: COLORS.combined, filter: (row) => row.ctype === "combined" && row.category === "All", value: (row) => row.y }
  ], { threshold: 0.5, label: "Overall mean target achievement trend" }, targetOption);
  renderCountryBars(document.querySelector("#features-chart"), data.countries, "nr_features");
  renderCountryBars(document.querySelector("#coverage-chart"), data.countries, "prop");
  renderHistogram(document.querySelector("#targets-chart"), data.targets.filter((row) => row.option === targetOption), targetOption);

  const areaTypes = ["natura2000", "cdda", "combined"];
  renderLineChart(document.querySelector("#area-chart"), data.areas, areaTypes.map((type) => ({
    label: type === "natura2000" ? "Natura 2000" : type === "cdda" ? "CDDA" : "Combined",
    color: COLORS[type], filter: (row) => row.ctype === type, value: (row) => row.area_km2
  })), { min: 0, max: 2000000, ticks: [0, 500000, 1000000, 1500000, 2000000], formatTick: (value) => `${(value / 1000000).toFixed(1)}m`, yAxisTitle: "Area (km²)", height: 300, label: "Cumulative protected area by area type" });

  const articles = [data.art12, data.art17Species, data.art17Habitats].flat().filter((row) => row.option === targetOption);
  const articleSeries = [
    { label: "Article 12", color: COLORS.art12, filter: (row) => row.art === "Art_12", value: (row) => row.y },
    { label: "Article 17", color: COLORS.art17, filter: (row) => row.art === "Art_17", value: (row) => row.y }
  ];
  [["species", "#article-species-chart"], ["habitats", "#article-habitats-chart"]].forEach(([category, selector]) => {
    const categoryRows = articles.filter((row) => String(row.category).toLowerCase() === category && row.ctype === "combined");
    renderTargetLineChart(document.querySelector(selector), categoryRows, articleSeries, { threshold: 0.5, label: `MTA by article for ${category}` }, targetOption);
  });
  ["species", "habitats"].forEach((category) => areaTypes.forEach((type) => {
    const selector = `#facet-${category}-${type}`;
    const categoryRows = articles.filter((row) => String(row.category).toLowerCase() === category && row.ctype === type);
    const availableSeries = articleSeries.filter((series) => categoryRows.some(series.filter));
    renderTargetLineChart(document.querySelector(selector), categoryRows, availableSeries, { threshold: 0.5, width: 520, height: 250, label: `MTA for ${category}, ${type}` }, targetOption);
  }));
}

/** Tracks the visible dashboard section and updates the active navigation link.
 * @returns {void} Registers an intersection observer for the navigation sections.
 */
function setupSectionNavigation() {
  const links = [...document.querySelectorAll(".section-nav a")];
  const sections = links.map((link) => document.querySelector(link.getAttribute("href"))).filter(Boolean);
  const observer = new IntersectionObserver((entries) => entries.forEach((entry) => {
    if (entry.isIntersecting) links.forEach((link) => link.classList.toggle("active", link.getAttribute("href") === `#${entry.target.id}`));
  }), { rootMargin: "-20% 0px -65% 0px" });
  sections.forEach((section) => observer.observe(section));
}

/** Adds per-chart expand controls and wires the shared inspection dialog.
 * @returns {void} Adds controls and dialog event listeners to the document.
 */
function setupChartExpansion() {
  const dialog = document.querySelector("#chart-dialog");
  const dialogTitle = document.querySelector("#chart-dialog-title");
  const dialogContent = document.querySelector("#chart-dialog-content");

  document.querySelectorAll(".chart-panel").forEach((panel, index) => {
    const heading = panel.querySelector(".panel-heading");
    if (!heading || heading.querySelector(".chart-expand")) return;
    const title = panel.querySelector("h2, h3, h4")?.textContent.trim() || `Chart ${index + 1}`;
    const button = document.createElement("button");
    button.className = "chart-expand";
    button.type = "button";
    button.textContent = "⛶";
    button.title = `Enlarge graph: ${title}`;
    button.setAttribute("aria-label", `Enlarge graph: ${title}`);
    button.setAttribute("aria-haspopup", "dialog");
    button.addEventListener("click", () => {
      const source = panel.querySelector(".chart-scroll");
      const svg = source?.querySelector("svg.chart-svg");
      if (!svg) return;

      dialogTitle.textContent = title;
      const enlargedChart = document.createElement("div");
      enlargedChart.className = "chart-dialog-chart";
      enlargedChart.append(svg.cloneNode(true));
      const legend = source.querySelector(".chart-legend");
      if (legend) enlargedChart.append(legend.cloneNode(true));
      dialogContent.replaceChildren(enlargedChart);
      dialog.showModal();
    });
    heading.append(button);
  });

  document.querySelector(".chart-dialog-close").addEventListener("click", () => dialog.close());
  dialog.addEventListener("click", (event) => {
    if (event.target === dialog) dialog.close();
  });
  updateChartExpansion();
}

/** Wires the further-limitations popup.
 * @returns {void} Adds open and close handlers to the limitations dialog.
 */
function setupLimitationsDialog() {
  const dialog = document.querySelector("#limitations-dialog");
  const openButton = document.querySelector("#limitations-open");
  const closeButton = document.querySelector(".limitations-dialog .chart-dialog-close");
  if (!dialog || !openButton || !closeButton) return;

  openButton.addEventListener("click", () => dialog.showModal());
  closeButton.addEventListener("click", () => dialog.close());
  dialog.addEventListener("click", (event) => {
    if (event.target === dialog) dialog.close();
  });
}

/** Disables expand controls for panels currently showing an empty state.
 * @returns {void} Synchronizes each button with its panel's SVG availability.
 */
function updateChartExpansion() {
  document.querySelectorAll(".chart-panel").forEach((panel) => {
    const button = panel.querySelector(".chart-expand");
    if (button) button.disabled = !panel.querySelector(".chart-scroll svg.chart-svg");
  });
}

/** Loads embedded data, renders the dashboard, and connects user controls.
 * @returns {void} Initializes dashboard rendering and event handling.
 */
function init() {
  try {
    document.querySelector("#data-version").textContent = window.MTA_VERSION || "unknown";
    if (!window.MTA_DATA) throw new Error("Dashboard data is missing. Rebuild dashboard-data.js from the export CSV files.");
    const targetSelect = document.querySelector("#target-select");
    renderDashboard(window.MTA_DATA, targetSelect.value);
    setupChartExpansion();
    setupLimitationsDialog();
    targetSelect.addEventListener("change", () => {
      renderDashboard(window.MTA_DATA, targetSelect.value);
      updateChartExpansion();
    });
    setupSectionNavigation();
  } catch (error) {
    const message = document.querySelector("#data-error");
    message.hidden = false;
    message.textContent = error.message;
    document.querySelectorAll(".chart-loading").forEach((loading) => { loading.textContent = "Data could not be loaded."; });
  }
}

init();