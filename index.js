// Import the ternary library
import TernaryLib from "./ternary.js";

// consts for the lattice view
const LATTICE_SVG_WIDTH = 500;
const LATTICE_SVG_HEIGHT = 500;
const ORIGIN_X = LATTICE_SVG_WIDTH / 2;
const ORIGIN_Y = LATTICE_SVG_HEIGHT / 2;
const SPACING_X = 35;
const SPACING_Y = -35;
const UNOCCUPIED_DOT_RADIUS = 7;
const EDGE_WIDTH = 2;

const GROUND_INDIGO = "#76f";

const statusElement = document.getElementById("status");

/**
 * Count occurrences of a character in a string
 */
function countChar(str, char) {
  return [...str].filter(c => c === char).length;
}
function displayStepVector(vector) {
  const keys = Object.keys(vector);
  if (keys.length === 0) return "0";
  const sizeIdentifiers = ["L", "m", "s"];
  return keys.map((k, i) => `${vector[k]}${sizeIdentifiers[i]}`).join("+");
}

// Use gcd and mod from TernaryLib
const { gcd, mod } = TernaryLib;

/**
 * Parse an equave ratio string (e.g., "2/1") and convert to cents
 * Returns 1200 (octave) if parsing fails
 */
function parseEquaveToCents(ratioStr) {
  const match = ratioStr.trim().match(/^(\d+)\/(\d+)$/);
  if (!match) {
    return 1200; // Default to octave
  }
  const numerator = parseInt(match[1], 10);
  const denominator = parseInt(match[2], 10);
  if (denominator === 0 || numerator <= 0 || denominator <= 0) {
    return 1200; // Default to octave
  }
  // cents = 1200 * log2(ratio)
  return 1200 * Math.log2(numerator / denominator);
}

/**
 * Get the current equave ratio string from the input field (normalized)
 * Returns { ratio: "m/n", num: m, den: n }
 */
function getEquaveRatio() {
  const input = document.getElementById("input-equave");
  if (!input) {
    return { ratio: "2/1", num: 2, den: 1 };
  }
  const value = input.value.trim();
  const match = value.match(/^(\d+)\/(\d+)$/);
  if (!match) {
    return { ratio: "2/1", num: 2, den: 1 };
  }
  const num = parseInt(match[1], 10);
  const den = parseInt(match[2], 10);
  const d = gcd(num, den);
  return { ratio: `${num / d}/${den / d}`, num: num / d, den: den / d };
}

/**
 * Get the current equave in cents from the input field
 */
function getEquaveCents() {
  const input = document.getElementById("input-equave");
  return input ? parseEquaveToCents(input.value) : 1200;
}

/**
 * Get the ED bound from the input field
 */
function getEdBound() {
  const input = document.getElementById("input-ed-bound");
  return input ? parseInt(input.value, 10) || TernaryLib.DEFAULT_ED_BOUND : TernaryLib.DEFAULT_ED_BOUND;
}

/**
 * Get the minimum s size in cents from the input field
 */
function getSLower() {
  const input = document.getElementById("input-s-lower");
  return input ? parseFloat(input.value) || TernaryLib.DEFAULT_S_LOWER_BOUND : TernaryLib.DEFAULT_S_LOWER_BOUND;
}

/**
 * Get the maximum s size in cents from the input field
 */
function getSUpper() {
  const input = document.getElementById("input-s-upper");
  return input ? parseFloat(input.value) || TernaryLib.DEFAULT_S_UPPER_BOUND : TernaryLib.DEFAULT_S_UPPER_BOUND;
}

function stepVectorLength(vector) {
  return Object.values(vector).reduce((sum, v) => sum + v, 0);
}

// Mutates `arr` a flat array matrix, swapping rows `i1` and `i2`.
function swapRows(arr, m, n, i1, i2) {
  if (i1 < 0 || i1 >= m || i2 < 0 || i2 >= m) {
    throw new Error(`swapRows(): matrix index out of bounds!`);
  }
  for (let j = 0; j < n; j++) {
    [arr[n * i1 + j], arr[n * i2 + j]] = [arr[n * i2 + j], arr[n * i1 + j]];
  }
}

// Scales row `i` of a matrix by `coeff`.
function multiplyRow(arr, m, n, i, coeff) {
  for (let j = 0; j < n; j++) {
    arr[n * i + j] *= coeff;
  }
}

// Does the operation "[row i2 of arr] += coeff * [row i1 of arr]".
function addMultipleOfFirstRowToSecond(arr, _, n, i1, i2, coeff) {
  for (let j = 0; j < n; j++) {
    arr[n * i2 + j] += coeff * arr[n * i1 + j];
  }
}

// Use Gaussian elimination to solve a linear system `left` * x = `right`.
// Both `left` and `right` are assumed to have `m` rows;
// `left` has `n1` columns, and `right` has `n2`.
// The function returns the RHS after this process.
function gaussianElimination(left, right, m, n1, n2) {
  // Don't mutate the input
  let leftClone = structuredClone(left);
  let rightClone = structuredClone(right);
  // Iterating over columns, clear all entries below the pivot
  for (let j = 0; j < Math.min(m, n1); ++j) {
    let i = j + 1;
    // Try to make the (j, j) entry nonzero by swapping rows
    while (i < m && Math.abs(leftClone[n1 * j + j]) < Number.EPSILON) {
      swapRows(leftClone, m, n1, j, i);
      swapRows(rightClone, m, n2, j, i);
      i++;
    }
    // Return null to indicate a singular matrix
    if (Math.abs(leftClone[n1 * j + j]) < Number.EPSILON) {
      return null;
    }
    // Clear lower triangle
    const pivot = leftClone[n1 * j + j];
    for (let i2 = j + 1; i2 < m; ++i2) {
      const target = leftClone[n1 * i2 + j];
      if (Math.abs(target) >= Number.EPSILON) {
        addMultipleOfFirstRowToSecond(leftClone, m, n1, j, i2, -target / pivot);
        addMultipleOfFirstRowToSecond(
          rightClone,
          m,
          n2,
          j,
          i2,
          -target / pivot,
        );
      }
    }
  }
  // Clear upper triangle
  for (let j = Math.min(m, n1) - 1; j >= 0; --j) {
    const pivot = leftClone[n1 * j + j];
    for (let i2 = 0; i2 < j; ++i2) {
      const target = leftClone[n1 * i2 + j];
      if (Math.abs(target) >= Number.EPSILON) {
        addMultipleOfFirstRowToSecond(leftClone, m, n1, j, i2, -target / pivot);
        addMultipleOfFirstRowToSecond(
          rightClone,
          m,
          n2,
          j,
          i2,
          -target / pivot,
        );
      }
    }
    // Scale rows so LHS gets 1 on diag
    // (Mutation can happen via another alias, though not via the `const` aliases we defined above.)
    multiplyRow(leftClone, m, n1, j, 1 / pivot);
    multiplyRow(rightClone, m, n2, j, 1 / pivot);
  }
  return rightClone;
}

// Makes a table in `tableElement` with the given `data`.
function makeTable(tableElement, data, header = "") {
  const tableViewTr = document.createElement("tr");
  let tableView = tableContent(data, header);
  tableViewTr.appendChild(tableView);
  tableElement.appendChild(tableViewTr);
}

// Return a new table view
function tableContent(data, header = "") {
  const table = tableHead(data, header);
  const tbody = table.createTBody();
  if (data[0] instanceof Array) {
    for (const [i] of data.entries()) {
      let row = tbody.insertRow();
      let cell1 = row.insertCell();
      cell1.appendChild(document.createTextNode(`${i + 1}`)); // row numbering
      for (const value of data[i].values()) {
        // iterate over columns
        let td = document.createElement("td");
        td.appendChild(document.createTextNode(value));
        row.appendChild(td);
      }
    }
  } else if (typeof data[0] === "object") {
    for (let i = 0; i < data.length; i++) {
      let row = tbody.insertRow();
      let cell1 = row.insertCell();
      cell1.appendChild(document.createTextNode(`${i + 1}`)); // row numbering
      for (const value of Object.values(data[i])) {
        // iterate over columns
        let td = document.createElement("td");
        td.appendChild(document.createTextNode(`${value}`));
        row.appendChild(td);
      }
    }
  } else {
    for (let i = 0; i < data.length; i++) {
      let row = tbody.insertRow();
      let cell1 = row.insertCell();
      cell1.appendChild(document.createTextNode(`${i + 1}`)); //row numbering
      let td = document.createElement("td");
      td.appendChild(document.createTextNode(data[i]));
      row.appendChild(td);
    }
  }
  return table;
}

function tableHead(data, header = "") {
  const table = document.createElement("table");
  const thead = table.createTHead();
  table.setAttribute("class", "scrollable clickable");
  let headRow = thead.insertRow();
  let th = document.createElement("th");
  th.appendChild(document.createTextNode("#"));
  headRow.appendChild(th);
  if (data) {
    if (data[0] instanceof Array) {
      for (let i = 0; i < data[0].length; ++i) {
        let th = document.createElement("th");
        th.appendChild(document.createTextNode(`${i}`));
        headRow.appendChild(th);
      }
    } else if (typeof data[0] === "object") {
      for (let key of Object.keys(data[0])) {
        let th = document.createElement("th");
        th.appendChild(document.createTextNode(`${key}`));
        headRow.appendChild(th);
      }
    } else {
      let th = document.createElement("th");
      th.appendChild(document.createTextNode(`${header}`));
      headRow.appendChild(th);
    }
  }
  return table;
}

// Main application code (using TernaryLib imported at top)
(function initApp() {
  // app state
  const appState = {
    word: null,
    latticeBasis: null,
    tuning: null,
    profile: null,
  };

  // New approach:
  // 1. draw a background 2D grid first
  // 2. represent x and y directions as generator and offset, whichever way fits better on the screen
  // 3. choose a zero point
  // TODO: Indicate what kind of guide frame it is. (simple, multiple/interleaved)
  //
  // TODO: show a legend for the different colored lines
  function createLatticeView(state, equave) {
    if (statusElement) {
      statusElement.innerText = "";

      let A;
      let B;
      let C;
      let D;

      // Create lattice visualization
      const latticeElement = document.getElementById("lattice-vis");
      if (latticeElement) {
        latticeElement.innerHTML = "";
        latticeElement.setAttribute("style", "vertical-align: text-top;");
        const svgTag = document.createElementNS(
          "http://www.w3.org/2000/svg",
          "svg",
        );
        svgTag.setAttribute("id", "lattice");
        svgTag.setAttribute("width", "100%");
        svgTag.setAttribute("height", "100%");
        const svgStyle = document.createElement("style");
        svgStyle.innerHTML = `
      .small {
        font: 20px sans-serif;
        fill: black;
      }`;
        if (state.word) {
          if (state.latticeBasis) {
            let g = state.latticeBasis[0];
            let h = state.latticeBasis[1];
            // Helper to get value from array or object with string keys
            const getVal = (v, i) =>
              Array.isArray(v) ? v[i] : (v[String(i)] ?? v[i] ?? 0);
            let sig = [0, 0, 0];
            let n = state.word.length;
            for (let i = 0; i < n; ++i) {
              switch (state.word[i]) {
                case "L":
                  ++sig[0];
                  break;
                case "m":
                  ++sig[1];
                  break;
                case "s":
                  ++sig[2];
                  break;
                default:
                  break;
              }
            }
            [A, B, C, D] = [1, 0, 0, 1];
            let [L_x, L_y, M_x, M_y, S_x, S_y] = gaussianElimination(
              [
                getVal(g, 0),
                getVal(g, 1),
                getVal(g, 2),
                getVal(h, 0),
                getVal(h, 1),
                getVal(h, 2),
                sig[0],
                sig[1],
                sig[2],
              ],
              [A, B, C, D, 0, 0],
              3,
              3,
              2,
            );
            for (let i = -20; i <= 20; ++i) {
              // Lines of constant g
              const p0x = ORIGIN_X + A * SPACING_X * i; // ith g offset
              const p0y = ORIGIN_X + B * SPACING_Y * i;
              const p1x = p0x + C * SPACING_X * 100; // add a large positive h offset
              const p1y = p0y + D * SPACING_Y * 100;
              const p2x = p0x + C * SPACING_X * -100; // add a large negative h offset
              const p2y = p0y + D * SPACING_Y * -100;
              // Lines of constant h
              const q0x = ORIGIN_X + C * SPACING_X * i;
              const q0y = ORIGIN_Y + D * SPACING_Y * i;
              const q1x = q0x + A * SPACING_X * 100;
              const q1y = q0y + B * SPACING_Y * 100;
              const q2x = q0x + A * SPACING_X * -100;
              const q2y = q0y + B * SPACING_Y * -100;

              // draw the line of constant g
              svgTag.innerHTML += `<line
                x1="${p1x}"
                y1="${p1y}"
                x2="${p2x}"
                y2="${p2y}"
                style="stroke:${GROUND_INDIGO}; stroke-width:${EDGE_WIDTH}"
              />`;
              // draw the line of constant h
              svgTag.innerHTML += `<line
                x1="${q1x}"
                y1="${q1y}"
                x2="${q2x}"
                y2="${q2y}"
                style="stroke:gray; stroke-width:${EDGE_WIDTH}"
              />`;
            }

            let currentX = ORIGIN_X;
            let currentY = ORIGIN_Y;
            // Track cumulative step counts for pitch calculation
            let stepCounts = { L: 0, m: 0, s: 0 };
            
            for (let deg = 0; deg < n; ++deg) {
              const { pitch, cents } = getPitchInfo(stepCounts, state.tuning, equave);
              const centsRounded = Math.round(cents);
              const tooltipText = `Degree ${deg}: ${pitch} (${centsRounded}¢)`;
              
              svgTag.innerHTML += `<g class="note-point" style="cursor: pointer;">
            <circle
              cx="${currentX}"
              cy="${currentY}"
              r="${UNOCCUPIED_DOT_RADIUS}"
              fill="white"
              stroke="white"
              stroke-width="1"
            >
              <title>${tooltipText}</title>
            </circle>
            <text
              x="${currentX - 3}"
              y="${currentY + 3}"
              fill="black"
              font-size="0.5em"
              style="pointer-events: none;"
            >${mod(deg, n)}</text>
          </g>`;
              switch (state.word[deg]) {
                case "L":
                  stepCounts.L++;
                  currentX += L_x * SPACING_X;
                  currentY += L_y * SPACING_Y; // SPACING_Y is negative since we represented y as positive.
                  break;
                case "m":
                  stepCounts.m++;
                  currentX += M_x * SPACING_X;
                  currentY += M_y * SPACING_Y;
                  break;
                case "s":
                  stepCounts.s++;
                  currentX += S_x * SPACING_X;
                  currentY += S_y * SPACING_Y;
                  break;
                default:
                  throw new Error("Invalid letter");
              }
            }
            // We deferred appending elements until now
            // Initial viewBox will be set by updateViewBox() below
            latticeElement.innerHTML += `<hr/><h2>Lattice view</h2><br/><small>Ternary scales are special in that they admit a JI-agnostic 2D lattice representation.<br/>Here the two dimensions g = ${alsoInCurrentTuning(g, state.tuning, equave)} and h = ${alsoInCurrentTuning(h, state.tuning, equave)} are two different generators. g is horizontal, h is vertical.</small>`;
            latticeElement.innerHTML += `<br/><small>Hover over the dots to see pitch information. Click and drag to pan, use mouse wheel or buttons to zoom.</small>`;
            // Add zoom buttons
            latticeElement.innerHTML += `<div class="controls">
        <button id="zoom-in">Zoom In (+)</button>
        <button id="zoom-out">Zoom Out (-)</button>
        <button id="reset-view">Reset View</button>
        <span style="margin-left: 20px;">Zoom: <span id="zoom-level">100%</span></span>
    </div>`;
            latticeElement.appendChild(svgTag);

            // Zoom functionality
            const zoomInButton = document.getElementById("zoom-in");
            const zoomOutButton = document.getElementById("zoom-out");
            const resetViewButton = document.getElementById("reset-view");
            const zoomLevelDisplay = document.getElementById("zoom-level");
            let viewBox = {
              x: 150,
              y: 150,
              width: LATTICE_SVG_WIDTH,
              height: LATTICE_SVG_HEIGHT,
            };

            let isPanning = false;
            let startPoint = { x: 0, y: 0 };
            let scale = 1.0;

            // Update viewBox attribute
            function updateViewBox() {
              svgTag.setAttribute(
                "viewBox",
                `${viewBox.x} ${viewBox.y} ${viewBox.width} ${viewBox.height}`,
              );
              zoomLevelDisplay.textContent = Math.round(scale * 100) + "%";
            }

            // Convert screen coordinates to SVG coordinates
            function getPointInSVG(e) {
              const CTM = svgTag.getScreenCTM();
              return {
                x: (e.clientX - CTM.e) / CTM.a,
                y: (e.clientY - CTM.f) / CTM.d,
              };
            }

            // Mouse down - start panning
            svgTag.addEventListener("mousedown", (e) => {
              isPanning = true;
              startPoint = getPointInSVG(e);
            });

            // Mouse move - pan
            svgTag.addEventListener("mousemove", (e) => {
              if (!isPanning) return;

              e.preventDefault();
              const currentPoint = getPointInSVG(e);
              const dx = currentPoint.x - startPoint.x;
              const dy = currentPoint.y - startPoint.y;

              viewBox.x -= dx;
              viewBox.y -= dy;

              updateViewBox();
            });

            // Mouse up - stop panning
            svgTag.addEventListener("mouseup", () => {
              isPanning = false;
            });

            svgTag.addEventListener("mouseleave", () => {
              isPanning = false;
            });

            // Wheel - zoom
            svgTag.addEventListener("wheel", (e) => {
              e.preventDefault();

              const point = getPointInSVG(e);
              const zoomFactor = e.deltaY < 0 ? 0.9 : 1.1;

              // Calculate new dimensions
              const newWidth = viewBox.width * zoomFactor;
              const newHeight = viewBox.height * zoomFactor;

              // Adjust position to zoom towards mouse cursor
              viewBox.x += (point.x - viewBox.x) * (1 - zoomFactor);
              viewBox.y += (point.y - viewBox.y) * (1 - zoomFactor);

              viewBox.width = newWidth;
              viewBox.height = newHeight;

              scale = 800 / viewBox.width;

              updateViewBox();
            });

            // Zoom by factor around center
            function zoomBy(factor) {
              const centerX = viewBox.x + viewBox.width / 2;
              const centerY = viewBox.y + viewBox.height / 2;
              viewBox.width *= factor;
              viewBox.height *= factor;
              viewBox.x = centerX - viewBox.width / 2;
              viewBox.y = centerY - viewBox.height / 2;
              scale = 800 / viewBox.width;
              updateViewBox();
            }

            // Button functions
            zoomInButton.addEventListener("click", () => zoomBy(0.8));
            zoomOutButton.addEventListener("click", () => zoomBy(1.25));
            resetViewButton.addEventListener("click", () => {
              viewBox = { x: 150, y: 150, width: LATTICE_SVG_WIDTH, height: LATTICE_SVG_HEIGHT };
              scale = 1;
              updateViewBox();
            });
            // Initialize
            updateViewBox();
          } else {
            // No lattice basis available - show a message instead of throwing
            latticeElement.innerHTML = `<hr/><h2>Lattice view</h2><br/><small>No suitable lattice basis found for this scale.</small>`;
          }
        } else {
          // No word selected
          latticeElement.innerHTML = "";
        }
      }
    }
  }

  // Function for showing the SonicWeave code
  function showSonicWeaveCode(state) {
    if (state.word) {
      const arity = new Set(Array.from(state.word)).size;
      if (state.tuning) {
        const element = document.getElementById("sw-code");
        if (element) {
          element.innerHTML = `<hr/><h2>SonicWeave code</h2>
        (for <a href="https://sw3.lumipakkanen.com/" target="_blank">Scale Workshop 3</a>)<br/>`;
          element.innerHTML += `<pre class="language-ocaml"><code class="language-ocaml" id="codeblock"></code></pre>`;

          // Make a "copy to clipboard" button
          const copyButtonLabel = "Copy";

          async function copyCode(block, button) {
            let code = block.querySelector("code");
            let text = code.innerText;

            await navigator.clipboard.writeText(text);

            // visual feedback that task is completed
            button.innerText = "Copied!";

            setTimeout(() => {
              button.innerText = copyButtonLabel;
            }, 700);
          }

          // use a class selector if available
          let blocks = document.querySelectorAll("pre");

          blocks.forEach((block) => {
            // only add button if browser supports Clipboard API
            if (navigator.clipboard) {
              let button = document.createElement("button");

              button.innerText = copyButtonLabel;
              block.appendChild(button);

              button.addEventListener("click", async () => {
                await copyCode(block, button);
              });
            }
          });

          let codeblock = document.getElementById("codeblock");
          const arr = Array.from(state.word);
          if (codeblock) {
            codeblock.innerHTML =
              arity === 3
                ? `let L = ${state.tuning[0]}
let m = ${state.tuning[1]}
let s = ${state.tuning[2]}
${arr.join(";")};
stack()`
                : arity === 2
                  ? `let L = ${state.tuning[0]}
let s = ${state.tuning[1]}
${arr.join(";")};
stack()`
                  : arity === 1
                    ? `let X = ${state.tuning[0]}
${arr.join(";")};
stack()`
                    : "Scales of rank > 3 are not supported";
          }
        }
      }
    }
  }

  function showScaleProfile(state, equave) {
    const el = document.getElementById("scale-profile");
    if (el) {
      el.innerHTML = "";
      const h2 = document.createElement("h2");
      h2.innerText = `Scale profile for ${state.word}`;
      el.appendChild(h2);
      if (state.profile) {
        // const ploidacot = state.profile["ploidacot"];
        const structure = state.profile["structure"];
        
        // Guide frame info (only if structure exists)
        if (structure) {
          el.innerHTML += `<b><a href="https://en.xen.wiki/w/Guide_frame
           " target="_blank">Guide frame</a></b><br/><small>`;
          let gsDisp =
            `${structure["gs"].map((g) => ` ${alsoInCurrentTuning(g, state.tuning, equave)}`)}`.slice(
              1,
            );
          el.innerHTML += `Guided <a href="https://en.xen.wiki/w/Generator_sequence" target="_blank">generator sequence</a> of ${stepVectorLength(structure["gs"][0])}-steps: GS(${gsDisp})<br/>`; // TODO prettify
          el.innerHTML += `Aggregate generator ${alsoInCurrentTuning(structure["aggregate"], state.tuning, equave)}<br/>`; // TODO prettify
          el.innerHTML += `Offsets ${structure["polyoffset"].map((g) => alsoInCurrentTuning(g, state.tuning, equave))}<br/>`; // TODO prettify
          el.innerHTML += `Multiplicity ${JSON.stringify(structure["multiplicity"])}<br/>`; // TODO prettify
          el.innerHTML += `Complexity ${JSON.stringify(structure["complexity"])}<br/><br/></small>`; // TODO prettify
        } else {
          el.innerHTML += `<b><a href="https://en.xen.wiki/w/Guide_frame" target="_blank">Guide frame</a></b><br/><small>No guide frame found.<br/><br/></small>`;
        }
        
        // Monotone MOS properties (always shown)
        el.innerHTML += `<b><a href="https://en.xen.wiki/w/Monotone-MOS_scale" target="_blank">Monotone MOS properties</a></b><br/><small>`;
        el.innerHTML += state.profile["lm"] ? `L = m<br/>` : "";
        el.innerHTML += state.profile["ms"] ? `m = s<br/>` : "";
        el.innerHTML += state.profile["s0"] ? `s = 0<br/>` : "";
        if (
          !state.profile["lm"] &&
          !state.profile["ms"] &&
          !state.profile["s0"]
        ) {
          el.innerHTML += `None<br/>`;
        }
        el.innerHTML += `<br/>`;
        
        // MOS substitution properties (always shown)
        const a = countChar(state.word, "L");
        const b = countChar(state.word, "m");
        const c = countChar(state.word, "s");

        el.innerHTML += `<b><a href="https://en.xen.wiki/w/MOS_substitution" target="_blank">MOS substitution</a> properties</b><br/>`;
        el.innerHTML += state.profile["subst_l_ms"]
          ? `subst ${a}L(${b}m${c}s)<br/>`
          : "";
        el.innerHTML += state.profile["subst_m_ls"]
          ? `subst ${b}m(${a}L${c}s)<br/>`
          : "";
        el.innerHTML += state.profile["subst_s_lm"]
          ? `subst ${c}s(${a}L${b}m)<br/>`
          : "";
        if (
          !state.profile["subst_l_ms"] &&
          !state.profile["subst_m_ls"] &&
          !state.profile["subst_s_lm"]
        ) {
          el.innerHTML += `None<br/>`;
        }

        // Chirality (always shown)
        if (state.profile["chirality"] === "Achiral") {
          el.innerHTML += `<br/><a href="https://en.xen.wiki/w/Chirality" target="_blank">Chirality</a>: Achiral`;
        } else {
          el.innerHTML += `<br/><a href="https://en.xen.wiki/w/Chirality" target="_blank">Chirality</a>: ${state.profile["chirality"]} (reversed: ${state.profile["reversed"]})`;
        }
        
        // Maximum variety (always shown)
        el.innerHTML += `<br/><br/><a href="https://en.xen.wiki/w/Maximum_variety" target="_blank">Maximum variety</a>: ${state.profile["mv"]}</small>`;
      }
    }
  }

  // display both the step vector in sum form and what interval it is in the current tuning
  function alsoInCurrentTuning(v, tuning, equave) {
    if (tuning["0"].includes("\\")) {
      if (equave.num !== 2 || equave.den !== 1) {
        const str0 = tuning["0"].substring(0, tuning["0"].indexOf("<"));
        const str1 = tuning["1"].substring(0, tuning["1"].indexOf("<"));
        const str2 = tuning["2"].substring(0, tuning["2"].indexOf("<"));
        const [numLstr, denLstr] = str0.split("\\");
        const [numMstr, _1] = str1.split("\\");
        const [numSstr, _2] = str2.split("\\");
        const numL = Number(numLstr);
        const numM = Number(numMstr);
        const numS = Number(numSstr);
        const ed = Number(denLstr);
        return `${displayStepVector(v)} (${numL * (v["0"] ?? 0) + numM * (v["1"] ?? 0) + numS * (v["2"] ?? 0)}\\${ed}<${equave.num}/${equave.den}>)`;
      } else {
        const [numLstr, denLstr] = tuning["0"].split("\\");
        const ed = Number(denLstr);
        const [numMstr, _3] = tuning["1"].split("\\");
        const [numSstr, _4] = tuning["2"].split("\\");
        const numL = Number(numLstr);
        const numM = Number(numMstr);
        const numS = Number(numSstr);
        return `${displayStepVector(v)} (${numL * (v["0"] ?? 0) + numM * (v["1"] ?? 0) + numS * (v["2"] ?? 0)}\\${ed})`;}
    } else if (tuning["0"].includes("/")) {
      const [numLstr, denLstr] = tuning["0"].split("/");
      const [numMstr, denMstr] = tuning["1"].split("/");
      const [numSstr, denSstr] = tuning["2"].split("/");
      const numL = Number(numLstr);
      const numM = Number(numMstr);
      const numS = Number(numSstr);
      const denL = Number(denLstr);
      const denM = Number(denMstr);
      const denS = Number(denSstr);
      const num =
        Math.pow(numL, v["0"] ?? 0) *
        Math.pow(numM, v["1"] ?? 0) *
        Math.pow(numS, v["2"] ?? 0);
      const den =
        Math.pow(denL, v["0"] ?? 0) *
        Math.pow(denM, v["1"] ?? 0) *
        Math.pow(denS, v["2"] ?? 0);
      const d = gcd(num, den);
      return `${displayStepVector(v)} (${num / d}/${den / d})`;
    } else {
      return displayStepVector(v);
    }
  }

  /**
   * Get pitch info for a scale degree given step counts
   * @param stepCounts - Object with counts of L, m, s steps: { L: n, m: n, s: n }
   * @param tuning - The current tuning object
   * @param equave - The equave { num, den, ratio }
   * @returns { pitch: string, cents: number }
   */
  function getPitchInfo(stepCounts, tuning, equave) {
    if (!tuning) return { pitch: "", cents: 0 };
    
    const nL = stepCounts.L || 0;
    const nM = stepCounts.m || 0;
    const nS = stepCounts.s || 0;
    
    // Calculate equave in cents
    const equaveCents = 1200 * Math.log2(equave.num / equave.den);
    
    if (tuning["0"].includes("\\")) {
      // ED tuning format like "5\12" or "5\12<3/1>"
      let str0 = tuning["0"];
      let str1 = tuning["1"];
      let str2 = tuning["2"];
      
      // Strip equave suffix if present
      if (str0.includes("<")) {
        str0 = str0.substring(0, str0.indexOf("<"));
        str1 = str1.substring(0, str1.indexOf("<"));
        str2 = str2.substring(0, str2.indexOf("<"));
      }
      
      const [numLstr, denLstr] = str0.split("\\");
      const [numMstr] = str1.split("\\");
      const [numSstr] = str2.split("\\");
      const stepsL = Number(numLstr);
      const stepsM = Number(numMstr);
      const stepsS = Number(numSstr);
      const ed = Number(denLstr);
      
      const totalSteps = nL * stepsL + nM * stepsM + nS * stepsS;
      const cents = (totalSteps / ed) * equaveCents;
      
      let pitch;
      if (equave.num !== 2 || equave.den !== 1) {
        pitch = `${totalSteps}\\${ed}<${equave.num}/${equave.den}>`;
      } else {
        pitch = `${totalSteps}\\${ed}`;
      }
      return { pitch, cents };
    } else if (tuning["0"].includes("/")) {
      // JI tuning format like "9/8"
      const [numLstr, denLstr] = tuning["0"].split("/");
      const [numMstr, denMstr] = tuning["1"].split("/");
      const [numSstr, denSstr] = tuning["2"].split("/");
      const numL = Number(numLstr);
      const numM = Number(numMstr);
      const numS = Number(numSstr);
      const denL = Number(denLstr);
      const denM = Number(denMstr);
      const denS = Number(denSstr);
      
      const num = Math.pow(numL, nL) * Math.pow(numM, nM) * Math.pow(numS, nS);
      const den = Math.pow(denL, nL) * Math.pow(denM, nM) * Math.pow(denS, nS);
      const d = gcd(num, den);
      const ratio = (num / d) / (den / d);
      const cents = 1200 * Math.log2(ratio);
      return { pitch: `${num / d}/${den / d}`, cents };
    } else {
      return { pitch: `${nL}L + ${nM}m + ${nS}s`, cents: 0 };
    }
  }

  function escapeHtml(text) {
    return text
      .replaceAll("&", "&amp;")
      .replaceAll("<", "&lt;")
      .replaceAll(">", "&gt;")
      .replaceAll(`'`, "&#39;")
      .replaceAll(`"`, "&quot;");
  }

  /**
   * Clear selection from both tuning tables and select a row
   */
  function selectTuningRow(jiTable, edTable, row) {
    jiTable.querySelector(`.selected`)?.classList.remove("selected");
    edTable.querySelector(`.selected`)?.classList.remove("selected");
    row.classList.add("selected");
  }

  // Helper to update all views
  function updateViews(equave) {
    showScaleProfile(appState, equave);
    createLatticeView(appState, equave);
    showSonicWeaveCode(appState);
  }

  const RANK_LIMIT = 3;
  const NO_HIGH_RANK_SCALES = `Scales of rank 0 or rank ${RANK_LIMIT + 1} or higher are not supported.`;

  const btnSig = document.getElementById("btn-sig");
  const btnWord = document.getElementById("btn-word");

  btnSig.addEventListener("click", () => {
    const sigQuery = document.getElementById("input-step-sig").value;
    const sig = `${sigQuery}`
      .split(" ")
      .map((str) => Number(str))
      .filter((n) => n);
    const arity = sig.filter((m) => m > 0).length;
    const scaleSize = sig.reduce((acc, m) => (acc += m), 0);
    try {
      if (scaleSize === 0) {
        statusElement.textContent = "Scale is empty. This is a bug.";
      } else if (arity > RANK_LIMIT) {
        statusElement.textContent = `Scales of rank 0 or rank ${RANK_LIMIT + 1} or higher are not supported.`;
      } else {
        if (arity > 0) {
          statusElement.textContent = "Computing...";

          document.getElementById("tables").innerHTML = `
      <hr /><h2>Tables</h2>
      <div
        style="
          overflow-y: auto;
          overflow-x: auto;
          vertical-align: top;
        "
      >
      <table>
        <tr>
                      <td>
                        Scales
                        <div
                          style="
                            overflow-y: auto;
                            overflow-x: auto;
                            vertical-align: top;
                            height: 420px;
                            width: 250px;
                          "
                        >
                          <table class="data" id="table-scales"></table>
                        </div>
                      </td>
                      <td>
                        JI tunings
                        <div
                          style="
                            overflow-y: auto;
                            overflow-x: auto;
                            vertical-align: top;
                            height: 420px;
                            width: 250px;
                          "
                        >
                          <table class="data" id="table-ji-tunings"></table>
                        </div>
                      </td>
                      <td>
                        ed(equave) tunings
                        <div
                          style="
                            overflow-y: auto;
                            overflow-x: auto;
                            vertical-align: top;
                            height: 420px;
                            width: 200px;
                          "
                        >
                          <table class="data" id="table-ed-tunings"></table>
                        </div>
                      </td></tr></table></div>`;
          const scaleTable = document.getElementById("table-scales");
          const jiTuningTable = document.getElementById("table-ji-tunings");
          const edTuningTable = document.getElementById("table-ed-tunings");
          const equave = getEquaveRatio();
          const sigResultData = TernaryLib.sigResult(sig, {
            lm: document.getElementById("monotone-lm").checked,
            ms: document.getElementById("monotone-ms").checked,
            s0: document.getElementById("monotone-s0").checked,
            ggsLen: Number(document.getElementById("ggs-len").value),
            ggsLenConstraint: document.querySelector(
              'input[name="ggs-len-constraint"]:checked',
            ).value,
            complexity: Number(document.getElementById("complexity").value),
            complexityConstraint: document.querySelector(
              'input[name="complexity-constraint"]:checked',
            ).value,
            mv: Number(document.getElementById("mv").value),
            mvConstraint: document.querySelector(
              'input[name="mv-constraint"]:checked',
            ).value,
            mosSubst: document.querySelector('input[name="mos-subst"]:checked')
              .value,
            equaveCents: getEquaveCents(),
            equaveRatio: equave.ratio,
            edBound: getEdBound(),
            sLower: getSLower(),
            sUpper: getSUpper(),
          });
          const scales = sigResultData["profiles"].map((j) => j["word"]);
          const latticeBases = sigResultData["profiles"].map(
            (j) => j["lattice_basis"],
          );
          const profiles = sigResultData["profiles"];

          const jiTunings = sigResultData["ji_tunings"];
          const edTunings = sigResultData["ed_tunings"];
          let letters;
          if (arity === 3) {
            letters = ["L", "m", "s"];
          } else if (arity === 2) {
            letters = ["L", "s"];
          } else if (arity === 1) {
            letters = ["X"];
          } else {
            letters = [...Array(arity).keys()].map((i) => `X${i}`);
          }
          statusElement.innerHTML = `<h1>Results for ${escapeHtml([...Array(arity).keys()].map((i) => `${sig[i]}${letters[i]}`).join(""))}</h1> (click on a table row to select a scale or a tuning)`;
          makeTable(scaleTable, scales, "scale");
          // add event listener for each non-head row
          const scaleRows = scaleTable.getElementsByTagName("tr");
          if (scaleRows.length >= 3) {
            scaleRows[2].classList.add("selected"); // For some reason 2 is the first row of a nonempty table.
            appState.word = scales[0];
            appState.profile = profiles[0];
            appState.latticeBasis = appState.profile["lattice_basis"];

            for (let i = 2; i < scaleRows.length; ++i) {
              scaleRows[i].addEventListener("click", async () => {
                // unselect selected row in either JI tuning table or ED tuning table
                scaleTable
                  .querySelector(`.selected`)
                  .classList.remove("selected");

                // select the row clicked on
                scaleRows[i].classList.add("selected");
                // get scale pattern
                let scaleWord = scales[i - 2];
                appState.profile = profiles[i - 2];
                appState.latticeBasis = latticeBases[i - 2];
                appState.word = scaleWord;
                updateViews(equave);
              });
            }
          }
          makeTable(jiTuningTable, jiTunings);
          const jiRows = jiTuningTable.getElementsByTagName("tr");
          for (let i = 2; i < jiRows.length; ++i) {
            const thisRow = jiRows[i];
            thisRow.addEventListener("click", () => {
              if (arity >= 1 && arity <= 3) {
                selectTuningRow(jiTuningTable, edTuningTable, thisRow);
                appState.tuning = jiTunings[i - 2];
                updateViews(equave);
              } else {
                statusElement.textContent = NO_HIGH_RANK_SCALES;
              }
            });
          }
          if (edTunings) {
            makeTable(edTuningTable, edTunings);
            const edRows = edTuningTable.getElementsByTagName("tr");
            edRows[2].classList.add("selected");
            const edTuning = edTunings[0];
            appState.tuning = edTuning;
            updateViews(equave);
            for (let i = 2; i < edRows.length; ++i) {
              const thisRow = edRows[i];
              thisRow.addEventListener("click", () => {
                selectTuningRow(jiTuningTable, edTuningTable, thisRow);
                const edTuning = edTunings[i - 2];
                appState.tuning = edTuning;
                updateViews(equave);
              });
            }
          }
          updateViews(equave);
        }
      }
    } catch (err) {
      statusElement.innerText = err;
    }
  });
  btnWord.addEventListener("click", () => {
    const query = document.getElementById("input-word").value;
    const arity = new Set(Array.from(query)).size;
    statusElement.textContent = "Computing...";
    const equave = getEquaveRatio();
    const wordResultData = TernaryLib.wordResult(query, {
      equaveCents: getEquaveCents(),
      equaveRatio: equave.ratio,
      edBound: getEdBound(),
      sLower: getSLower(),
      sUpper: getSUpper(),
    });
    const profile = wordResultData["profile"];
    const brightestMode = wordResultData["profile"]["word"];
    const jiTunings = wordResultData["ji_tunings"];
    const edTunings = wordResultData["ed_tunings"];
    document.getElementById("tables").innerHTML = `
                <td>
                  JI tunings
                  <div
                    style="
                      overflow-y: auto;
                      overflow-x: auto;
                      vertical-align: top;
                      height: 420px;
                      width: 250px;
                    "
                  >
                    <table class="data" id="table-ji-tunings"></table>
                  </div>
                </td>
                <td>
                  ed(equave) tunings
                  <div
                    style="
                      overflow-y: auto;
                      overflow-x: auto;
                      vertical-align: top;
                      height: 420px;
                      width: 200px;
                    "
                  >
                    <table class="data" id="table-ed-tunings"></table>
                  </div>
                </td>`;
    const jiTuningTable = document.getElementById("table-ji-tunings");
    const edTuningTable = document.getElementById("table-ed-tunings");
    appState.word = brightestMode;
    try {
      statusElement.innerHTML = `<h1>Results for ${appState.word}</h1>(click on a table row to select a tuning)`;

      makeTable(jiTuningTable, jiTunings);
      const jiRows = jiTuningTable.getElementsByTagName("tr");
      for (let i = 2; i < jiRows.length; ++i) {
        const thisRow = jiRows[i];
        thisRow.addEventListener("click", () => {
          if (arity >= 1 && arity <= 3) {
            selectTuningRow(jiTuningTable, edTuningTable, thisRow);
            appState.tuning = jiTunings[i - 2];
            updateViews(equave);
          } else {
            statusElement.textContent = NO_HIGH_RANK_SCALES;
          }
        });
      }
      if (edTunings) {
        makeTable(edTuningTable, edTunings);
        const edRows = edTuningTable.getElementsByTagName("tr");
        edRows[2].classList.add("selected");
        for (let i = 2; i < edRows.length; ++i) {
          const thisRow = edRows[i];
          thisRow.addEventListener("click", () => {
            selectTuningRow(jiTuningTable, edTuningTable, thisRow);
            const edTuning = edTunings[i - 2];
            appState.tuning = edTuning;
            updateViews(equave);
          });
        }
      }
      appState.tuning = edTunings[0];
      appState.profile = profile;
      appState.latticeBasis = appState.profile["lattice_basis"];
      updateViews(equave);
    } catch (err) {
      statusElement.innerText = err;
    }
  });
})();
