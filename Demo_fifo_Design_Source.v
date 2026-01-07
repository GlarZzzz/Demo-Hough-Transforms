`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 01/07/2026 03:10:23 PM
// Design Name: 
// Module Name: Demo_fifo_Design_Source
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module Demo_fifo_Design_Source #(

    parameter DATA_WIDTH = 21,   
    parameter DEPTH      = 65536/4,
    
    parameter ADDR_WIDTH = $clog2(DEPTH)
    )
    
    (
    
    input  wire                  clk,
    input  wire                  resetn, 
    input  wire [DATA_WIDTH-1:0] wr_data,
    input  wire                  wr_en,
    output wire                  full,
    output wire [DATA_WIDTH-1:0] rd_data,
    input  wire                  rd_en,
    output wire                  empty
    
    );
    
        reg [DATA_WIDTH-1:0] mem_array[0:DEPTH-1];

   
    reg [ADDR_WIDTH-1:0] wr_ptr;
    reg [ADDR_WIDTH-1:0] rd_ptr;
    
    reg [ADDR_WIDTH:0]   item_count; 

    
    reg [DATA_WIDTH-1:0] rd_data_reg;

    assign rd_data = rd_data_reg;
    assign full    = (item_count == DEPTH);
    assign empty   = (item_count == 0);
    
    
    always @(posedge clk) begin
        if (!resetn) begin 
        end else begin
            if (wr_en && !full) begin
                mem_array[wr_ptr] <= wr_data;
            end
        end
    end

    
    always @(posedge clk) begin
        if (!resetn) begin 
            wr_ptr      <= 0;
            rd_ptr      <= 0;
            item_count  <= 0;
            rd_data_reg <= 0; 
        end else begin
            
            if (wr_en && !full) begin
                wr_ptr <= wr_ptr + 1;
            end

            
            if (rd_en && !empty) begin
                
                rd_ptr <= rd_ptr + 1;
               
                rd_data_reg <= mem_array[rd_ptr]; 
            end

            
            case ({wr_en && !full, rd_en && !empty})
                2'b01: 
                    item_count <= item_count - 1;
                2'b10: 
                    item_count <= item_count + 1;
                2'b11: 
                    item_count <= item_count; 
                default: 
                    item_count <= item_count;
            endcase
        end
    end
    
endmodule
