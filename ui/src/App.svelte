<script lang="ts">
  import "./app.css";
  import Navbar from "$lib/components/Navbar.svelte";
  import { Label } from "$lib/components/ui/label";
  import { Separator } from "$lib/components/ui/separator/index.js";
  import { Input } from "$lib/components/ui/input/index.js";
  import * as Card from "$lib/components/ui/card/index.js";
  import * as Select from "$lib/components/ui/select/index.js";
  import * as HoverCard from "$lib/components/ui/hover-card/index.js";
  import { cn } from "$lib/utils"; // shadcn helper for conditional classes
  import FileUpload from "$lib/components/FileUpload.svelte";

  import logo from "/src/assets/logo.svg";
  import Progress from "$lib/components/ui/progress/progress.svelte";

  let peptide_file: File[];
  let prog_value = 20;
  let currentStep = 1;

  const steps = [
    { id: 1, label: "Upload Peptide Network" },
    { id: 2, label: "Omic Data Selection" },
    { id: 3, label: "Analysis" },
    { id: 4, label: "Results" },
  ];

  async function check_network() {
    const name = peptide_file.name;
    document.getElementById("status").textContent = "Starting task...";
    document.getElementById("result").textContent = "";

    const response = await fetch("/start_task", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name }),
    });

    const data = await response.json();
    currentTask = data.task_id;

    document.getElementById("status").textContent =
      "Task started. Checking progress...";
    pollTask();
  }
</script>

<Navbar />

<logo><img src={logo} alt="Logo" class="mx-auto mb-4 w-32" /></logo>

<div class="container mx-auto p-4 prose prose-lg max-w-3/4">
  <h1 class="text-center">XLRanker</h1>
  <div class="grid grid-cols-[200px_1fr] gap-6">
    <!-- Sidebar -->
    <aside class="flex flex-col space-y-1 text-sm">
      {#each steps as step}
        <div
          class={cn(
            "rounded-lg px-3 py-1 transition-colors",
            step.id === currentStep && "font-bold text-black",
            step.id < currentStep && "text-gray-700",
            step.id > currentStep && "text-gray-400"
          )}
        >
          {step.label}
        </div>
      {/each}
    </aside>
    <div>
      <div
        class="w-full grid space-x-4 text-sm border px-4 rounded-lg bg-white shadow-md"
      >
        <div class="px-4">
          <h3>Peptide Network</h3>
          Please upload a <HoverCard.Root>
            <HoverCard.Trigger>Peptide Network</HoverCard.Trigger>
            <HoverCard.Content>
              <div class="grid gap-2">
                <h4 class="font-semibold">Example</h4>
                <pre id="example" class="whitespace-pre-wrap">
EQLAKLMATLR	IIGLDQVAGM
EQLAKLMATLR	KGMFR
EQLAKLMATLR	KGMFR
EQLDNQLDAK	SLNLKHIK
EQLDNQLDAY	SLNLKHIK
</pre>
              </div>
            </HoverCard.Content>
          </HoverCard.Root> file containing peptide pairs, one pair per line. Each
          line should contain two peptides separated by a tab or comma character.

          <div class="grid gap-2">
            <h3>Upload Network</h3>
            <FileUpload bind:files={peptide_file} />
            <h4>Peptide Separator</h4>
            <Select.Root type="single">
              <Select.Trigger>Select separator</Select.Trigger>
              <Select.Content>
                <Select.Item value="tab">Tab</Select.Item>
                <Select.Item value="comma">Comma</Select.Item>
                <Select.Item value="space">Space</Select.Item>
              </Select.Content>
            </Select.Root>
          </div>

          <h3>Sequence Mapping</h3>
          <p>
            To map peptides to their corresponding proteins, please upload a
            FASTA file or use the default.
          </p>
        </div>
      </div>
      <Progress value={prog_value} max={100} id="progress" class="my-4" />
    </div>
  </div>
</div>
