// The Sudoku page: fill the grid, prove in the browser, verify, tamper.
const PUZZLE = "530070000600195000098000060800060003400803001700020006060000280000419005000080079";
const SOLUTION = "534678912672195348198342567859761423426853791713924856961537284287419635345286179";

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const grid = $<HTMLDivElement>("grid");
const status = $<HTMLParagraphElement>("status");
const proofBox = $<HTMLElement>("proof");
const worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });

const cells: HTMLInputElement[] = [];
for (let i = 0; i < 81; i++) {
  const inp = document.createElement("input");
  inp.inputMode = "numeric";
  inp.maxLength = 1;
  inp.setAttribute("aria-label", `row ${Math.floor(i / 9) + 1}, column ${(i % 9) + 1}`);
  if (PUZZLE[i] !== "0") {
    inp.value = PUZZLE[i];
    inp.readOnly = true;
    inp.className = "given";
  }
  inp.addEventListener("input", () => {
    inp.value = inp.value.replace(/[^1-9]/g, "");
    proof = null;
    render();
  });
  if (Math.floor(i / 9) % 3 === 2 && i < 72) inp.classList.add("bottom");
  if (i % 9 % 3 === 2 && i % 9 !== 8) inp.classList.add("right");
  grid.appendChild(inp);
  cells.push(inp);
}

let proof: Uint8Array | null = null;
let pub: string[] = [];

function render() {
  $<HTMLButtonElement>("verify").disabled = !proof;
  $<HTMLButtonElement>("tamper").disabled = !proof;
  proofBox.textContent = proof ? Array.from(proof, (b) => b.toString(16).padStart(2, "0")).join("") : "";
  $("proof-size").textContent = proof ? `${proof.length} bytes` : "";
}

function say(text: string, kind: "" | "ok" | "bad" = "") {
  status.textContent = text;
  status.className = kind;
}

$("fill").onclick = () => {
  cells.forEach((c, i) => (c.value = SOLUTION[i]));
  proof = null;
  render();
};

$("prove").onclick = () => {
  const sol = cells.map((c) => c.value || "0");
  const inputs = `puzzle: ${PUZZLE.split("").join(" ")}\nsolution: ${sol.join(" ")}\n`;
  say("Proving in WebAssembly…");
  worker.postMessage({ type: "prove", inputs });
};

$("verify").onclick = () => {
  if (!proof) return;
  say("Verifying in TypeScript…");
  worker.postMessage({ type: "verify", proof, pub });
};

$("tamper").onclick = () => {
  if (!proof) return;
  const i = Math.floor(Math.random() * proof.length);
  proof = proof.slice();
  proof[i] ^= 1 << Math.floor(Math.random() * 8);
  render();
  say(`Flipped one bit of byte ${i}. Verify again.`);
};

worker.onmessage = (ev: MessageEvent) => {
  const m = ev.data;
  if (m.type === "proved") {
    if (m.ok) {
      proof = m.proof;
      pub = m.pub;
      say(`Proof made in ${(m.ms / 1000).toFixed(1)} s. It says the grid is solved, without saying how.`, "ok");
    } else {
      proof = null;
      say(m.error.startsWith("Unsatisfied") ? "That is not a solution of this puzzle: there is nothing to prove." : `Error: ${m.error}`, "bad");
    }
    render();
  } else if (m.type === "verified") {
    say(m.ok ? `Valid: whoever made this proof knows a solution (checked in ${m.ms.toFixed(0)} ms).` : "Invalid proof.", m.ok ? "ok" : "bad");
  } else {
    say(`Error: ${m.error}`, "bad");
  }
};
render();
