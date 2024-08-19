// consts for the lattice view
const LATTICE_SVG_WIDTH = 600;
const LATTICE_SVG_HEIGHT = 600;
const ORIGIN_X = LATTICE_SVG_WIDTH / 2;
const ORIGIN_Y = LATTICE_SVG_HEIGHT / 2;
const SPACING_X = 35;
const SPACING_Y = -35;
const UNOCCUPIED_DOT_RADIUS = 7;
const EDGE_WIDTH = 2;

const GROUND_INDIGO = "#6c00da";

// global state (TODO: replace this with arguments for functions that show scale and tuning data)
let currentWord = null;
let currentLatticeBasis = null;
let currentTuning = null;
let currentProfile = null;

const statusElement = document.getElementById("status");
function displayStepVector(vector) {
  const keys = [...Object.keys(vector)];
  const sizeIdentifiers = ["L", "M", "s"];
  if (keys.length === 0) {
    return "0";
  }
  let result = "";
  for (let i = 0; i < keys.length; ++i) {
    result += `${vector[keys[i]]}${sizeIdentifiers[i]}`;
    if (i < keys.length - 1) {
      result += "+";
    }
  }
  return result;
}

function gcd(a, b) {
  if (a === 0) {
    return b;
  } else if (b === 0) {
    return a;
  } else {
    return a > b ? gcd(a % b, b) : a === b ? a : gcd(a, b % a);
  }
}

// Modulo that works like Python %
function modulo(a, m) {
  return (m + (a % m)) % m;
}

function stepVectorLength(vector) {
  let result = 0;
  for (let key of Object.keys(vector)) {
    result += vector[key];
  }
  return result;
}
// The m x m identity matrix.
function eye(m) {
  let arr = Array(m * m);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < m; j++) {
      arr[m * i + j] = i === j ? 1 : 0;
    }
  }
  return arr;
}

// Mutates `arr` a flat array matrix, swapping rows `i1` and `i2`.
function swapRows(arr, m, n, i1, i2) {
  if (i1 < 0 || i1 >= m) {
    throw new Error(
      `switchRows(): matrix index out of bounds! i1: ${i1} but m: ${m}`,
    );
  }
  if (i2 < 0 || i2 >= m) {
    throw new Error(
      `switchRows(): matrix index out of bounds! i2: ${i2} but m: ${m}`,
    );
  }
  for (let j = 0; j < n; j++) {
    const new1 = arr[n * i2 + j];
    const new2 = arr[n * i1 + j];
    [arr[n * i1 + j], arr[n * i2 + j]] = [new1, new2];
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
    // Clear lower triangle{
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

// Finds the inverse of `matrix`,
// If `matrix` is not invertible, throws an error.
function inverseMatrix(matrix, m) {
  return gaussianElimination(matrix, eye(m), m, m, m);
}

// Makes a table in `tableElement` with the given `data`.
function makeTable(tableElement, data, header = "") {
  const tableViewTr = document.createElement("tr");
  let tableView = tableContent(data, header);
  tableViewTr.appendChild(tableView);
  tableElement.appendChild(tableViewTr);
}

// Return a new table view
// TODO: make table entries selectable
function tableContent(data, header = "") {
  const table = tableHead(data, header);
  const tbody = table.createTBody();
  if (data[0] instanceof Array) {
    for (const [i, element] of data.entries()) {
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
    for (const [i, element] of data.entries()) {
      let row = tbody.insertRow();
      let cell1 = row.insertCell();
      cell1.appendChild(document.createTextNode(`${i + 1}`)); // row numbering
      for (const value of Object.values(element)) {
        // iterate over columns
        let td = document.createElement("td");
        td.appendChild(document.createTextNode(`${value}`));
        row.appendChild(td);
      }
    }
  } else {
    for (const [i, _] of data.entries()) {
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

import("./pkg").then((wasm) => {
  // New approach:
  // 1. draw a background 2D grid first
  // 2. represent x and y directions as generator and offset, whichever way fits better on the screen
  // 3. choose a zero point
  // Indicate what kind of guide frame it is. (simple, interleaved, multiple)
  //
  // TODO: show a legend for the different colored lines
  function createLatticeView() {
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
        const svgStyle = document.createElement("style");
        svgStyle.innerHTML = `
      .small {
        font: 20px sans-serif;
        fill: black;
      }`;
        if (currentWord) {
          if (currentLatticeBasis) {
            let g = currentLatticeBasis[0];
            let h = currentLatticeBasis[1];
            let sig = [0, 0, 0];
            let n = currentWord.length;
            for (let i = 0; i < n; ++i) {
              switch (currentWord[i]) {
                case "L":
                  ++sig[0];
                  break;
                case "M":
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
              [g[0], g[1], g[2], h[0], h[1], h[2], sig[0], sig[1], sig[2]],
              [A, B, C, D, 0, 0],
              3,
              3,
              2,
            );
            for (let i = -8; i <= 8; ++i) {
              // Lines of constant g
              const p0x = ORIGIN_X + A * SPACING_X * i; // ith g offset
              const p0y = ORIGIN_X + B * SPACING_Y * i;
              const p1x = p0x + C * SPACING_X * 10; // add a large positive h offset
              const p1y = p0y + D * SPACING_Y * 10;
              const p2x = p0x + C * SPACING_X * -10; // add a large negative h offset
              const p2y = p0y + D * SPACING_Y * -10;
              // Lines of constant h
              const q0x = ORIGIN_X + C * SPACING_X * i;
              const q0y = ORIGIN_Y + D * SPACING_Y * i;
              const q1x = q0x + A * SPACING_X * 10;
              const q1y = q0y + B * SPACING_Y * 10;
              const q2x = q0x + A * SPACING_X * -10;
              const q2y = q0y + B * SPACING_Y * -10;

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
            for (let deg = 0; deg < n; ++deg) {
              svgTag.innerHTML += `<circle 
            cx="${currentX}"
            cy="${currentY}"
            r="${UNOCCUPIED_DOT_RADIUS}"
            color="white"
            stroke="black" 
            stroke-width="1"
          />
          <text
            x="${currentX - 3}"
            y="${currentY + 3}"
            fill="white"
            font-size="0.5em"
          >${modulo(deg, n)}</text>`;
              switch (currentWord[deg]) {
                case "L":
                  currentX += L_x * SPACING_X;
                  currentY += L_y * SPACING_Y; // SPACING_Y is negative since we represented y as positive.
                  break;
                case "M":
                  currentX += M_x * SPACING_X;
                  currentY += M_y * SPACING_Y;
                  break;
                case "s":
                  currentX += S_x * SPACING_X;
                  currentY += S_y * SPACING_Y;
                  break;
                default:
                  throw new Error("Invalid letter");
              }
            }
            // We deferred appending elements until now
            svgTag.setAttribute(
              "viewBox",
              `0 0 ${LATTICE_SVG_WIDTH} ${LATTICE_SVG_HEIGHT}`,
            );
            latticeElement.innerHTML += `<h2>Lattice view</h2><br/><small>Ternary scales are special in that they admit a JI-agnostic 2D lattice representation.<br/>Here the two dimensions g = ${alsoInCurrentTuning(g)} and h = ${alsoInCurrentTuning(h)} are two different generators. g is horizontal, h is vertical.</small>`;
            latticeElement.appendChild(svgTag);
          } else {
            throw new Error("No suitable lattice basis");
          }
        } else {
          throw new Error("`currentWord` is null or undefined");
        }
      }
    }
  }

  // Function for showing the SonicWeave code
  function showSonicWeaveCode() {
    if (currentWord) {
      const arity = new Set(Array.from(currentWord)).size;
      if (currentTuning) {
        const element = document.getElementById("sw-code");
        if (element) {
          element.innerHTML = `<h2>SonicWeave code</h2>
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
          const arr = Array.from(currentWord);
          if (codeblock) {
            codeblock.innerHTML =
              arity === 3
                ? `let L = ${currentTuning[0]}
let M = ${currentTuning[1]}
let s = ${currentTuning[2]}
${arr.join(";")};
stack()`
                : arity === 2
                  ? `let L = ${currentTuning[0]}
let s = ${currentTuning[1]}
${arr.join(";")};
stack()`
                  : arity === 1
                    ? `let X = ${currentTuning[0]}
${arr.join(";")};
stack()`
                    : "Scales of rank > 3 are not supported";
          }
        }
      }
    }
  }

  function showScaleProfile() {
    const el = document.getElementById("scale-profile");
    if (el) {
      el.innerHTML = "";
      const h2 = document.createElement("h2");
      h2.innerText = `Scale profile for ${currentWord}`;
      el.appendChild(h2);
      if (currentProfile) {
        const structure = currentProfile["structure"];
        if (structure) {
          el.innerHTML += `<b>Guide frame</b><br/><small>`;
          let gsDisp =
            `${structure["gs"].map((g) => ` ${alsoInCurrentTuning(g)}`)}`.slice(
              1,
            );
          el.innerHTML += `Guided generator sequence of ${stepVectorLength(structure["gs"][0])}-steps: GS(${gsDisp})<br/>`; // TODO prettify
          el.innerHTML += `Aggregate generator ${alsoInCurrentTuning(structure["aggregate"])}<br/>`; // TODO prettify
          el.innerHTML += `Interleaving polyoffset ${structure["polyoffset"].map((g) => alsoInCurrentTuning(g))}<br/>`; // TODO prettify
          el.innerHTML += `Multiplicity ${JSON.stringify(structure["multiplicity"])}<br/>`; // TODO prettify
          el.innerHTML += `Complexity ${JSON.stringify(structure["complexity"])}<br/><br/>`; // TODO prettify
          el.innerHTML += `<b>Monotone MOS properties</b><br/>`;
          el.innerHTML += currentProfile["lm"] ? `L = M<br/>` : "";
          el.innerHTML += currentProfile["ms"] ? `M = s<br/>` : "";
          el.innerHTML += currentProfile["s0"] ? `s = 0<br/>` : "";
          if (
            !currentProfile["lm"] &&
            !currentProfile["ms"] &&
            !currentProfile["s0"]
          ) {
            el.innerHTML += `None<br/>`;
          }
          el.innerHTML += `<br/>Maximum variety ${currentProfile["mv"]}</small>`;
        }
      }
    }
  }

  // display both the step vector in sum form and what interval it is in the current tuning
  function alsoInCurrentTuning(v) {
    if (currentTuning["0"].includes("/")) {
      const [numLstr, denLstr] = currentTuning["0"].split("/");
      const [numMstr, denMstr] = currentTuning["1"].split("/");
      const [numSstr, denSstr] = currentTuning["2"].split("/");
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
    } else if (currentTuning["0"].includes("\\")) {
      const [numLstr, denLstr] = currentTuning["0"].split("\\");
      const ed = Number(denLstr);
      const [numMstr, _1] = currentTuning["1"].split("\\");
      const [numSstr, _2] = currentTuning["2"].split("\\");
      const numL = Number(numLstr);
      const numM = Number(numMstr);
      const numS = Number(numSstr);
      return `${displayStepVector(v)} (${numL * (v["0"] ?? 0) + numM * (v["1"] ?? 0) + numS * (v["2"] ?? 0)}\\${ed})`;
    } else {
      return displayStepVector(v);
    }
  }

  function escapeHtml(text) {
    text.replaceAll("&", "&amp;");
    text.replaceAll("<", "&lt;");
    text.replaceAll(">", "&gt;");
    text.replaceAll(`'`, "&#39;");
    text.replaceAll(`"`, "&quot;");
    return text;
  }
  const RANK_LIMIT = 3;
  const NO_HIGH_RANK_SCALES = `Scales of rank 0 or rank ${RANK_LIMIT + 1} or higher are not supported.`;

  btnSig = document.getElementById("btn-sig");
  btnWord = document.getElementById("btn-word");

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
                        edo tunings
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
          const scaleTable = document.getElementById("table-scales");
          const jiTuningTable = document.getElementById("table-ji-tunings");
          const edoTuningTable = document.getElementById("table-ed-tunings");
          const sigResult = wasm.sig_result(
            sig,
            document.getElementById("monotone-lm").checked,
            document.getElementById("monotone-ms").checked,
            document.getElementById("monotone-s0").checked,
            Number(document.getElementById("ggs-len").value),
            document.querySelector('input[name="ggs-len-constraint"]:checked')
              .value,
            Number(document.getElementById("complexity").value),
            document.querySelector(
              'input[name="complexity-constraint"]:checked',
            ).value,
            Number(document.getElementById("mv").value),
            document.querySelector('input[name="mv-constraint"]:checked').value,
            document.querySelector('input[name="mos-subst"]:checked').value,
          );
          const scales = sigResult["profiles"].map((j) => j["word"]);
          const latticeBases = sigResult["profiles"].map(
            (j) => j["lattice_basis"],
          );
          console.log(latticeBases);

          const profiles = sigResult["profiles"];

          const jiTunings = sigResult["ji_tunings"];
          const edoTunings = sigResult["ed_tunings"];
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
            currentWord = scales[0];
            currentProfile = profiles[0];
            currentLatticeBasis = currentProfile["lattice_basis"];

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
                currentProfile = profiles[i - 2];
                currentLatticeBasis = latticeBases[i - 2];
                currentWord = scaleWord;
                showSonicWeaveCode();
                showScaleProfile();

                try {
                  createLatticeView();
                } catch (err) {
                  statusElement.innerText = err;
                }
              });
            }
          }
          makeTable(jiTuningTable, jiTunings);
          const jiRows = jiTuningTable.getElementsByTagName("tr");
          for (let i = 2; i < jiRows.length; ++i) {
            const thisRow = jiRows[i];
            thisRow.addEventListener("click", () => {
              if (arity >= 1 && arity <= 3) {
                // unselect selected row in either JI tuning table or ED tuning table
                if (jiTuningTable.querySelector(`.selected`)) {
                  jiTuningTable
                    .querySelector(`.selected`)
                    .classList.remove("selected");
                }
                if (edoTuningTable.querySelector(`.selected`)) {
                  edoTuningTable
                    .querySelector(`.selected`)
                    .classList.remove("selected");
                }
                // select the row clicked on
                thisRow.classList.add("selected");
                currentTuning = jiTunings[i - 2];
                showSonicWeaveCode();
                showScaleProfile();
                createLatticeView();
              } else {
                statusElement.textContent = NO_HIGH_RANK_SCALES;
              }
            });
          }
          if (edoTunings) {
            makeTable(edoTuningTable, edoTunings);
            const edRows = edoTuningTable.getElementsByTagName("tr");
            edRows[2].classList.add("selected");
            const edoTuning = edoTunings[0];
            currentTuning = edoTuning;
            showScaleProfile();
            createLatticeView(currentWord, currentProfile["structure"]);
            showSonicWeaveCode();
            for (let i = 2; i < edRows.length; ++i) {
              const thisRow = edRows[i];
              thisRow.addEventListener("click", () => {
                // unselect selected row in either JI tuning table or ED tuning table
                if (jiTuningTable.querySelector(`.selected`)) {
                  jiTuningTable
                    .querySelector(`.selected`)
                    .classList.remove("selected");
                }
                if (edoTuningTable.querySelector(`.selected`)) {
                  edoTuningTable
                    .querySelector(`.selected`)
                    .classList.remove("selected");
                }
                // select the row clicked on
                thisRow.classList.add("selected");
                const edoTuning = edoTunings[i - 2];
                currentTuning = edoTuning;
                showScaleProfile();
                createLatticeView();
                showSonicWeaveCode();
              });
            }
          }
          showSonicWeaveCode();
          showScaleProfile();
          createLatticeView();
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
    const wordResult = wasm.word_result(query);
    const profile = wordResult["profile"];
    const brightestMode = wordResult["profile"]["word"];
    const jiTunings = wordResult["ji_tunings"];
    const edoTunings = wordResult["ed_tunings"];
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
                  edo tunings
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
    const edoTuningTable = document.getElementById("table-ed-tunings");
    currentWord = brightestMode;
    try {
      statusElement.innerHTML = `<h1>Results for ${currentWord}</h1>(click on a table row to select a tuning)`;

      const sig_ = { L: 0, M: 0, s: 0 };
      for (let i = 0; i < currentWord.length; ++i) {
        sig_[currentWord[i]] += 1;
      }
      const sig = Object.values(sig_); // values will be extracted ordered by the keys; for indexing by 0, 1, 2... rather than L, m, s

      makeTable(jiTuningTable, jiTunings);
      const jiRows = jiTuningTable.getElementsByTagName("tr");
      for (let i = 2; i < jiRows.length; ++i) {
        const thisRow = jiRows[i];
        thisRow.addEventListener("click", () => {
          if (arity >= 1 && arity <= 3) {
            // unselect selected row in either JI tuning table or ED tuning table
            if (jiTuningTable.querySelector(`.selected`)) {
              jiTuningTable
                .querySelector(`.selected`)
                .classList.remove("selected");
            }
            if (edoTuningTable.querySelector(`.selected`)) {
              edoTuningTable
                .querySelector(`.selected`)
                .classList.remove("selected");
            }
            // select the row clicked on
            thisRow.classList.add("selected");
            currentTuning = jiTunings[i - 2];
            showScaleProfile();
            createLatticeView();
            showSonicWeaveCode();
          } else {
            statusElement.textContent = NO_HIGH_RANK_SCALES;
          }
        });
      }
      if (edoTunings) {
        makeTable(edoTuningTable, edoTunings);
        const edRows = edoTuningTable.getElementsByTagName("tr");
        edRows[2].classList.add("selected");
        for (let i = 2; i < edRows.length; ++i) {
          const thisRow = edRows[i];
          thisRow.addEventListener("click", () => {
            // unselect selected row in either JI tuning table or ED tuning table
            if (jiTuningTable.querySelector(`.selected`)) {
              jiTuningTable
                .querySelector(`.selected`)
                .classList.remove("selected");
            }
            if (edoTuningTable.querySelector(`.selected`)) {
              edoTuningTable
                .querySelector(`.selected`)
                .classList.remove("selected");
            }
            // select the row clicked on
            thisRow.classList.add("selected");
            const edoTuning = edoTunings[i - 2];
            currentTuning = edoTuning;
            showScaleProfile();
            createLatticeView();
            showSonicWeaveCode();
          });
        }
      }
      const edoTuning = edoTunings[0];
      const ed = edoTuning["ed"];
      currentTuning = edoTuning;
      showSonicWeaveCode();
      currentProfile = profile;
      currentLatticeBasis = currentProfile["lattice_basis"];
      const guideFrame = currentProfile["structure"];
      showScaleProfile();
      createLatticeView();
      showSonicWeaveCode();
    } catch (err) {
      statusElement.innerText = err;
    }
  });
});
