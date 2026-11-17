// A static file server for the demo: node serve.mjs, then http://localhost:8000
import { createServer } from "node:http";
import { readFile } from "node:fs/promises";
import { extname, join, normalize } from "node:path";

const root = import.meta.dirname;
const types = { ".html": "text/html; charset=utf-8", ".js": "text/javascript", ".wasm": "application/wasm", ".bin": "application/octet-stream", ".circuit": "text/plain" };
createServer(async (req, res) => {
  const path = normalize(decodeURIComponent(new URL(req.url, "http://x").pathname)).replace(/^[\\/]+/, "");
  const file = join(root, path || "index.html");
  if (!file.startsWith(root)) return res.writeHead(403).end();
  try {
    const body = await readFile(file);
    res.writeHead(200, { "content-type": types[extname(file)] ?? "application/octet-stream" }).end(body);
  } catch {
    res.writeHead(404).end("not found");
  }
}).listen(Number(process.env.PORT ?? 8000), () => console.log(`http://localhost:${process.env.PORT ?? 8000}`));
